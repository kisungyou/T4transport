#' Distance between Histograms
#' 
#' Compute the \eqn{p}-Wasserstein distance between two 1D histograms that share 
#' the same binning, i.e., same breaks. The histograms are treated as discrete 
#' probability measures supported at bin midpoints with masses given by 
#' normalized counts. Uses the exact 1D monotone OT algorithm, not LP nor entropic regularization.
#' 
#' @param hist1 a histogram object (class \code{"histogram"}).
#' @param hist2 a histogram object (class \code{"histogram"}) with the same breaks as \code{hist1}.
#' @param p an exponent for the order of the distance (default: 2).
#' 
#' @return a named list containing\describe{
#' \item{distance}{\eqn{\mathcal{W}_p} distance value.}
#' }
#' 
#' @examples 
#' \donttest{
#' #----------------------------------------------------------------------
#' #                      Binned from Gaussian and Uniform
#' #
#' # Create two types of histograms with the same binning. One is from 
#' # the standard normal and the other from uniform distribution in [-5,5].
#' #----------------------------------------------------------------------
#' # GENERATE 20 HISTOGRAMS
#' set.seed(100)
#' hist20 = list()
#' bk = seq(from=-10, to=10, length.out=20) # common breaks
#' for (i in 1:10){
#'   hist20[[i]] = hist(stats::rnorm(100), breaks=bk, plot=FALSE)
#'   hist20[[i+10]] = hist(stats::runif(100, min=-5, max=5), breaks=bk, plot=FALSE)
#' }
#' 
#' # COMPUTE THE PAIRWISE DISTANCE
#' pdmat = array(0,c(20,20))
#' for (i in 1:19){
#'   for (j in (i+1):20){
#'     pdmat[i,j] = histdist(hist20[[i]], hist20[[j]], p=2)$distance
#'     pdmat[j,i] = pdmat[i,j]
#'   }
#' }
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(pty="s")
#' image(pdmat, axes=FALSE, main="Pairwise 2-Wasserstein Distance between Histograms")
#' par(opar)
#' }
#' 
#' @concept histogram
#' @export
histdist <- function(hist1, hist2, p=2){
  # INPUTS
  if (!inherits(hist1, "histogram")){
    stop("* histdist: input 'hist1' should be a histogram object.")
  }
  if (!inherits(hist2, "histogram")){
    stop("* histdist: input 'hist2' should be a histogram object.")
  }
  if (is.null(hist1$breaks) || is.null(hist2$breaks)) {
    stop("* histdist: both inputs must have a 'breaks' component (histogram objects).")
  }
  if (is.null(hist1$counts) || is.null(hist2$counts)) {
    stop("* histdist: both inputs must have a 'counts' component.")
  }
  
  # check equal breaks
  breaks1 <- hist1$breaks
  breaks2 <- hist2$breaks
  
  # enforce same binning
  if (length(breaks1) != length(breaks2) || any(breaks1 != breaks2)) {
    stop("* histdist: 'hist1' and 'hist2' must have identical 'breaks'.")
  }
  
  ## bin midpoints
  mids <- if (!is.null(hist1$mids)) {
    hist1$mids
  } else {
    (breaks1[-1L] + breaks1[-length(breaks1)]) / 2
  }
  
  counts1 <- as.numeric(hist1$counts)
  counts2 <- as.numeric(hist2$counts)
  
  if (any(counts1 < 0) || any(counts2 < 0)) {
    stop("* histdist: negative counts are not allowed.")
  }
  
  sum1 <- sum(counts1)
  sum2 <- sum(counts2)
  if (sum1 <= 0 || sum2 <= 0) {
    stop("* histdist: both histograms must have positive total mass.")
  }
  
  ## normalize to probability masses
  w1 <- counts1 / sum1
  w2 <- counts2 / sum2
  
  ## order p >= 1
  p <- max(1, as.numeric(p))
  K <- length(mids)
  
  ## 1D optimal transport by greedy mass matching
  i <- 1L
  j <- 1L
  mass1 <- w1[i]
  mass2 <- w2[j]
  cost  <- 0.0
  
  tol <- sqrt(.Machine$double.eps)
  
  while (i <= K && j <= K) {
    m <- min(mass1, mass2)
    if (m > 0) {
      d <- abs(mids[i] - mids[j])
      cost <- cost + m * (d ^ p)
    }
    
    mass1 <- mass1 - m
    mass2 <- mass2 - m
    
    if (mass1 <= tol) {
      i <- i + 1L
      if (i <= K) {
        mass1 <- w1[i]
      }
    }
    if (mass2 <= tol) {
      j <- j + 1L
      if (j <= K) {
        mass2 <- w2[j]
      }
    }
  }
  
  output = list()
  output$distance = cost^(1 / p)
  return(output)
}