#' Barycenter of Histograms 
#'
#' @description
#' Given multiple histograms represented as \code{"histogram"} S3 objects, compute
#' their 2-Wasserstein barycenter using the exact 1D quantile characterization. 
#' All input histograms must have identical breaks. 
#'
#' @param hists a length-\eqn{N} list of histograms (\code{"histogram"} objects)
#'   of same breaks.
#' @param weights a weight for each histogram; if \code{NULL} (default), uniform
#'   weights are used. Otherwise, it should be a length-\eqn{N} vector of
#'   nonnegative weights.
#' @param L number of quantile levels used to approximate the barycenter
#'   (default: 2000). Larger \code{L} gives a more accurate approximation at
#'   increased computational cost.
#'
#' @return a \code{"histogram"} object representing the Wasserstein barycenter.
#'   
#' @examples 
#' \donttest{
#' #----------------------------------------------------------------------
#' #                      Binned from Two Gaussians
#' #
#' # EXAMPLE : Very Small Example for CRAN; just showing how to use it!
#' #----------------------------------------------------------------------
#' # GENERATE FROM TWO GAUSSIANS WITH DIFFERENT MEANS
#' set.seed(100)
#' x  = stats::rnorm(1000, mean=-4, sd=0.5)
#' y  = stats::rnorm(1000, mean=+4, sd=0.5)
#' bk = seq(from=-10, to=10, length.out=20)
#' 
#' # HISTOGRAMS WITH COMMON BREAKS
#' histxy = list()
#' histxy[[1]] = hist(x, breaks=bk, plot=FALSE)
#' histxy[[2]] = hist(y, breaks=bk, plot=FALSE)
#' 
#' # COMPUTE
#' hh = histbary(histxy)
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' barplot(histxy[[1]]$density, col=rgb(0,0,1,1/4), 
#'         ylim=c(0, 0.75), main="Two Histograms")
#' barplot(histxy[[2]]$density, col=rgb(1,0,0,1/4), 
#'         ylim=c(0, 0.75), add=TRUE)
#' barplot(hh$density, main="Barycenter",
#'         ylim=c(0, 0.75))
#' par(opar)
#' }
#'
#' @concept histogram
#' @export
histbary <- function(hists, weights = NULL, L = 2000L) {
  name.f <- "histbary"
  check.f <- check_hists(hists, name.f)
  # check_hists is assumed to:
  #   - validate 'hists'
  #   - ensure same breaks
  #   - provide midpoints in check.f$midpts
  #   - (possibly) provide densities in check.f$density (not used here)
  
  nhists <- length(hists)
  if (nhists < 1L) {
    stop("* histbary: 'hists' must be a non-empty list of histogram objects.")
  }
  
  # reference histogram
  href <- hists[[1L]]
  if (!inherits(href, "histogram") ||
      is.null(href$breaks) || is.null(href$counts)) {
    stop("* histbary: each element of 'hists' must be a 'histogram' with 'breaks' and 'counts'.")
  }
  
  base_breaks <- href$breaks
  K <- length(base_breaks) - 1L
  p = 2.0
  
  #-------------------------------------------------
  # Barycenter weights
  #-------------------------------------------------
  myweights <- valid_multiple_weight(weights, nhists, name.f)
  myweights <- myweights / base::sum(myweights)
  
  #-------------------------------------------------
  # Bin midpoints and CDFs of each histogram
  #-------------------------------------------------
  mids <- check.f$midpts
  if (length(mids) != K) {
    stop("* histbary: internal error - length of midpoints and breaks do not match.")
  }
  
  CDFs <- matrix(0, nrow = nhists, ncol = K)
  
  for (i in seq_len(nhists)) {
    hi <- hists[[i]]
    ci <- as.numeric(hi$counts)
    
    if (any(ci < 0)) {
      stop("* histbary: negative counts are not allowed.")
    }
    total_i <- sum(ci)
    if (total_i <= 0) {
      stop("* histbary: each histogram must have positive total mass.")
    }
    
    pi <- ci / total_i                 # probability mass per bin
    CDFs[i, ] <- cumsum(pi)           # CDF at right edge of each bin
  }
  
  #-------------------------------------------------
  # Approximate barycenter quantile function
  #-------------------------------------------------
  L <- as.integer(L)
  if (L <= 0L) {
    stop("* histbary: 'L' must be a positive integer.")
  }
  
  bary_samples <- numeric(L)
  
  for (ell in seq_len(L)) {
    t <- (ell - 0.5) / L  # quantile level in (0,1)
    
    x_vals <- numeric(nhists)
    for (i in seq_len(nhists)) {
      j <- which(CDFs[i, ] >= t)[1L]
      if (is.na(j)) j <- K
      x_vals[i] <- mids[j]
    }
    
    # quantile of barycenter at level t: weighted average of input quantiles
    bary_samples[ell] <- sum(myweights * x_vals)
  }
  
  #-------------------------------------------------
  # Re-bin barycenter samples onto original breaks
  #-------------------------------------------------
  hb_raw <- graphics::hist(bary_samples, breaks = base_breaks, plot = FALSE)
  
  # approximate bin probabilities
  prob_hat <- hb_raw$counts / sum(hb_raw$counts)
  
  # match total count scale of the first histogram
  total_ref <- sum(href$counts)
  new_counts <- round(total_ref * prob_hat)
  
  # densities so that integral = 1
  bin_widths <- diff(base_breaks)
  new_density <- prob_hat / bin_widths
  
  #-------------------------------------------------
  # Wrap and return as a histogram object
  #-------------------------------------------------
  output <- href
  output$counts  <- new_counts
  output$density <- new_density
  output$xname   <- "barycenter"
  class(output)  <- "histogram"
  
  return(output)
}