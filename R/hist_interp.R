#' Interpolation between Histograms 
#'
#' Given two histograms represented as \code{"histogram"} S3 objects with
#' identical breaks, compute interpolated histograms along the 2-Wasserstein
#' geodesic connecting them. In 1D, this is achieved by linear interpolation
#' of quantile functions (displacement interpolation).
#'
#' @param hist1 a histogram (\code{"histogram"} object).
#' @param hist2 another histogram with the same \code{breaks} as \code{hist1}.
#' @param t a scalar or numeric vector in \eqn{[0,1]} specifying interpolation
#'   times. \code{t = 0} returns \code{hist1}, \code{t = 1} returns \code{hist2}.
#' @param L number of quantile levels used to approximate the geodesic
#'   (default: 2000). Larger \code{L} gives a more accurate approximation at
#'   increased computational cost.
#'
#' @return
#' If \code{length(t) == 1}, a single \code{"histogram"} object representing the
#' interpolated distribution at time \code{t}.
#' If \code{length(t) > 1}, a length-\code{length(t)} list of \code{"histogram"}
#' objects.
#'
#' @examples
#' \donttest{
#' #----------------------------------------------------------------------
#' #                      Interpolating Two Gaussians
#' #
#' # The source histogram is created from N(-5,1/4).
#' # The target histogram is created from N(+5,4)
#' #----------------------------------------------------------------------
#' # SETTING
#' set.seed(123)
#' x_source = rnorm(1000, mean=-5, sd=1/2)
#' x_target = rnorm(1000, mean=+5, sd=2)
#' 
#' # BUILD HISTOGRAMS WITH COMMON BREAKS
#' bk = seq(from=-8, to=12, by=2)
#' h1 = hist(x_source, breaks=bk, plot=FALSE)
#' h2 = hist(x_target, breaks=bk, plot=FALSE)
#' 
#' # INTERPOLATE WITH 5 GRID POINTS
#' h_path <- histinterp(h1, h2, t = seq(0, 1, length.out = 8))
#' 
#' # VISUALIZE
#' y_slim <- c(0, max(h1$density, h2$density)) # shared y-limit
#' xt     <- round(h1$mids, 1) # x-ticks
#' 
#' opar <- par(no.readonly = TRUE)
#' par(mfrow = c(2,4), pty = "s")
#' for (i in 1:8){
#'   if (i < 2){
#'     barplot(h_path[[i]]$density,  names.arg=xt, ylim=y_slim,
#'             main="Source", col=rgb(0,0,1,1/4))
#'   } else if (i > 7){
#'     barplot(h_path[[i]]$density,  names.arg=xt, ylim=y_slim,
#'             main="Target", col=rgb(1,0,0,1/4))
#'   } else {
#'     barplot(h_path[[i]]$density,  names.arg=xt, ylim=y_slim, 
#'             col="gray90", main=sprintf("t = %.3f", (i-1)/7))
#'   }
#' }
#' par(opar)
#' }
#'
#' @concept histogram
#' @export
histinterp <- function(hist1, hist2, t = 0.5, L = 2000L) {
  name.f <- "histinterp"
  h1 <- hist1
  h2 <- hist2
  
  #-------------------------------------------------
  # 0. Basic checks via existing helper
  #-------------------------------------------------
  hlist  <- list(h1, h2)
  check.f <- check_hists(hlist, name.f)
  # check_hists is assumed to:
  #   - validate 'hists'
  #   - ensure same breaks
  #   - provide midpoints in check.f$midpts
  
  # reference histogram
  href <- h1
  if (!inherits(href, "histogram") ||
      is.null(href$breaks) || is.null(href$counts)) {
    stop("* histinterp: 'hist1' must be a 'histogram' with 'breaks' and 'counts'.")
  }
  
  base_breaks <- href$breaks
  K           <- length(base_breaks) - 1L
  
  #-------------------------------------------------
  # 1. t in [0,1]
  #-------------------------------------------------
  t <- as.numeric(t)
  if (any(t < 0 | t > 1)) {
    stop("* histinterp: all 't' values must lie in [0,1].")
  }
  
  #-------------------------------------------------
  # 2. Bin midpoints and CDFs of each histogram
  #-------------------------------------------------
  mids <- check.f$midpts
  if (length(mids) != K) {
    stop("* histinterp: internal error - length of midpoints and breaks do not match.")
  }
  
  # build CDFs for h1 and h2
  build_cdf <- function(h) {
    cts <- as.numeric(h$counts)
    if (any(cts < 0)) {
      stop("* histinterp: negative counts are not allowed.")
    }
    tot <- sum(cts)
    if (tot <= 0) {
      stop("* histinterp: each histogram must have positive total mass.")
    }
    p <- cts / tot
    cumsum(p)
  }
  
  C1 <- build_cdf(h1)
  C2 <- build_cdf(h2)
  
  #-------------------------------------------------
  # 3. Quantile vectors for h1 and h2
  #-------------------------------------------------
  L <- as.integer(L)
  if (L <= 0L) {
    stop("* histinterp: 'L' must be a positive integer.")
  }
  
  q1 <- numeric(L)
  q2 <- numeric(L)
  
  for (ell in seq_len(L)) {
    u <- (ell - 0.5) / L  # quantile level
    
    j1 <- which(C1 >= u)[1L]
    if (is.na(j1)) j1 <- K
    q1[ell] <- mids[j1]
    
    j2 <- which(C2 >= u)[1L]
    if (is.na(j2)) j2 <- K
    q2[ell] <- mids[j2]
  }
  
  #-------------------------------------------------
  # 4. For each t, interpolate quantiles and re-bin
  #-------------------------------------------------
  make_hist_at_t <- function(tau) {
    # Displacement interpolation in quantile space
    q_tau <- (1 - tau) * q1 + tau * q2
    
    # Treat q_tau as L samples from the interpolated distribution
    h_raw <- graphics::hist(q_tau, breaks = base_breaks, plot = FALSE)
    
    # approximate bin probabilities
    prob_hat <- h_raw$counts / sum(h_raw$counts)
    
    # match total count scale of the first histogram
    total_ref <- sum(href$counts)
    new_counts <- round(total_ref * prob_hat)
    
    # densities so that integral = 1
    bin_widths <- diff(base_breaks)
    new_density <- prob_hat / bin_widths
    
    out <- href
    out$counts  <- new_counts
    out$density <- new_density
    out$xname   <- sprintf("interpolation t=%.3f", tau)
    class(out)  <- "histogram"
    out
  }
  
  if (length(t) == 1L) {
    return(make_hist_at_t(t))
  } else {
    out_list <- vector("list", length(t))
    for (k in seq_along(t)) {
      out_list[[k]] <- make_hist_at_t(t[k])
    }
    return(out_list)
  }
}