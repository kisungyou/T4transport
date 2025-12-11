#' Wasserstein Median of Histograms
#'
#' Given multiple histograms represented as \code{"histogram"} S3 objects with
#' common breaks, compute their Fréchet (geometric) median under the
#' 2-Wasserstein distance. In 1D, this is implemented by mapping histograms
#' to their quantile functions and running a Weiszfeld-type algorithm for 
#' the geometric median in the Hilbert space of quantile
#' functions.
#'
#' @param hists a length-\eqn{N} list of histograms (\code{"histogram"} objects)
#'   with identical \code{breaks}.
#' @param weights a weight for each histogram; if \code{NULL} (default), uniform
#'   weights are used. Otherwise, it should be a length-\eqn{N} vector of
#'   nonnegative weights.
#' @param L number of quantile levels used to approximate the median
#'   (default: 2000). Larger \code{L} gives a more accurate approximation at
#'   increased computational cost.
#' @param ... extra parameters including \describe{
#'   \item{abstol}{stopping criterion for iterations (default: 1e-8).}
#'   \item{maxiter}{maximum number of iterations (default: 496).}
#'   \item{print.progress}{logical; whether to show current iteration
#'        (default: \code{FALSE}).}
#' }
#'
#' @return a \code{"histogram"} object representing the Wasserstein median.
#' 
#' @examples 
#' \donttest{
#' #----------------------------------------------------------------------
#' #                      Binned from Two Gaussians
#' #
#' # Generate 12 histograms from N(-4,1/4) and 8 from N(4,1/4)
#' #----------------------------------------------------------------------
#' # COMMON SETTING
#' set.seed(100)
#' bk = seq(from=-10, to=10, length.out=20)
#' n_signal = 12
#' 
#' # GENERATE HISTOGRAMS WITH COMMON BREAKS
#' hist_all = list()
#' for (i in 1:n_signal){
#'   hist_all[[i]] = hist(stats::rnorm(200, mean=-4, sd=0.5), breaks=bk)
#' }
#' for (j in (n_signal+1):20){
#'   hist_all[[j]] = hist(stats::rnorm(200, mean=+4, sd=0.5), breaks=bk)
#' }
#' 
#' # COMPUTE THE BARYCENTER AND THE MEDIAN 
#' h_bary = histbary(hist_all)
#' h_med  = histmed(hist_all)
#' 
#' # VISUALIZE
#' xt   <- round(h_med$mids, 1) 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' barplot(hist_all[[1]]$density, col=rgb(0,0,1,1/4), 
#'         ylim=c(0, 0.75), main="Two Types", names.arg=xt)
#' barplot(hist_all[[20]]$density, col=rgb(1,0,0,1/4), 
#'         ylim=c(0, 0.75), add=TRUE)
#' barplot(h_med$density,  names.arg=xt, main="Median", ylim=c(0, 0.75))
#' barplot(h_bary$density, names.arg=xt, main="Barycenter", ylim=c(0, 0.75))
#' par(opar)
#' }
#' 
#' @concept histogram
#' @export
histmed <- function(hists, weights = NULL, L = 2000L, ...) {
  name.f  <- "histmed"
  
  #-------------------------------------------------
  # 0. Basic checks via existing helper
  #-------------------------------------------------
  check.f <- check_hists(hists, name.f)
  # check_hists is assumed to:
  #   - validate 'hists'
  #   - ensure same breaks
  #   - provide midpoints in check.f$midpts
  
  nhists <- length(hists)
  if (nhists < 1L) {
    stop("* histmed: 'hists' must be a non-empty list of histogram objects.")
  }
  
  href <- hists[[1L]]
  if (!inherits(href, "histogram") ||
      is.null(href$breaks) || is.null(href$counts)) {
    stop("* histmed: each element of 'hists' must be a 'histogram' with 'breaks' and 'counts'.")
  }
  
  base_breaks <- href$breaks
  K          <- length(base_breaks) - 1L
  
  #-------------------------------------------------
  # 1. Enforce p = 2 (W2 metric)
  #-------------------------------------------------
  p <- as.double(2.0)
  if (abs(p - 2) > sqrt(.Machine$double.eps)) {
    stop("* histmed: W2-median is currently implemented only for p = 2.")
  }
  
  #-------------------------------------------------
  # 2. Weights
  #-------------------------------------------------
  myweights <- valid_multiple_weight(weights, nhists, name.f)
  myweights <- myweights / base::sum(myweights)
  
  #-------------------------------------------------
  # 3. Bin midpoints and CDFs of each histogram
  #-------------------------------------------------
  mids <- check.f$midpts
  if (length(mids) != K) {
    stop("* histmed: internal error - length of midpoints and breaks do not match.")
  }
  
  CDFs <- matrix(0, nrow = nhists, ncol = K)
  
  for (i in seq_len(nhists)) {
    hi <- hists[[i]]
    ci <- as.numeric(hi$counts)
    
    if (any(ci < 0)) {
      stop("* histmed: negative counts are not allowed.")
    }
    total_i <- sum(ci)
    if (total_i <= 0) {
      stop("* histmed: each histogram must have positive total mass.")
    }
    
    pi <- ci / total_i                 # probability mass per bin
    CDFs[i, ] <- cumsum(pi)           # CDF at right edge of each bin
  }
  
  #-------------------------------------------------
  # 4. Build quantile vectors Y[i,ell] ≈ F_i^{-1}(t_ell)
  #-------------------------------------------------
  L <- as.integer(L)
  if (L <= 0L) {
    stop("* histmed: 'L' must be a positive integer.")
  }
  
  Y <- matrix(0, nrow = nhists, ncol = L)
  
  for (ell in seq_len(L)) {
    t <- (ell - 0.5) / L  # quantile level in (0,1)
    for (i in seq_len(nhists)) {
      j <- which(CDFs[i, ] >= t)[1L]
      if (is.na(j)) j <- K
      Y[i, ell] <- mids[j]
    }
  }
  
  #-------------------------------------------------
  # 5. Weiszfeld iteration for geometric median in R^L
  #-------------------------------------------------
  params <- list(...)
  pnames <- names(params)
  
  myiter <- max(
    1L,
    round(ifelse(("maxiter" %in% pnames), params$maxiter, 496))
  )
  mytol  <- max(
    100 * .Machine$double.eps,
    as.double(ifelse(("abstol" %in% pnames), params$abstol, 1e-8))
  )
  myshow <- as.logical(ifelse(("print.progress" %in% pnames), params$print.progress, FALSE))
  
  # initial guess: weighted mean of quantile vectors (i.e., barycenter quantile)
  z_cur <- as.numeric(myweights %*% Y)  # length L
  
  for (it in seq_len(myiter)) {
    # distances from current estimate
    diff_mat <- Y - matrix(z_cur, nrow = nhists, ncol = L, byrow = TRUE)
    dists    <- sqrt(rowSums(diff_mat^2))
    
    # handle potential zero distances (exact coincidence)
    # if any exactly zero, that point is already a minimizer in theory
    if (any(dists < mytol)) {
      # pick the coinciding point as median
      idx_zero <- which(dists < mytol)[1L]
      z_cur <- Y[idx_zero, ]
      if (myshow) {
        message(sprintf("* histmed: converged exactly at sample %d.\n", idx_zero))
      }
      break
    }
    
    # Weiszfeld weights: w_i / ||z - y_i||
    w_eff <- myweights / dists
    denom <- sum(w_eff)
    if (denom <= 0) {
      # degenerate case; fall back to previous
      warning("* histmed: Weiszfeld denominator nonpositive; returning last iterate.")
      break
    }
    
    z_new <- as.numeric((w_eff %*% Y) / denom)
    
    # check convergence in L2 norm
    step_norm <- sqrt(sum((z_new - z_cur)^2))
    if (myshow) {
      message(sprintf("* histmed: iter %d, step norm = %.3e", it, step_norm))
    }
    z_cur <- z_new
    
    if (step_norm < mytol) {
      if (myshow) {
        message(sprintf("* histmed: converged in %d iterations.", it))
      }
      break
    }
  }
  
  z_med <- z_cur  # final quantile vector
  
  #-------------------------------------------------
  # 6. Re-bin median quantiles onto original breaks
  #-------------------------------------------------
  # We can treat z_med as L samples from the median distribution
  med_samples <- z_med
  
  hmed_raw <- graphics::hist(med_samples, breaks = base_breaks, plot = FALSE)
  
  # approximate bin probabilities
  prob_hat <- hmed_raw$counts / sum(hmed_raw$counts)
  
  # match total count scale of the first histogram
  total_ref <- sum(href$counts)
  new_counts <- round(total_ref * prob_hat)
  
  # densities so that integral = 1
  bin_widths <- diff(base_breaks)
  new_density <- prob_hat / bin_widths
  
  #-------------------------------------------------
  # 7. Wrap and return as a histogram object
  #-------------------------------------------------
  output <- href
  output$counts  <- new_counts
  output$density <- new_density
  output$xname   <- "median"
  class(output)  <- "histogram"
  return(output)
}