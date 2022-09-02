#' Barycenter of Histograms by Benamou et al. (2015)
#' 
#' Given multiple histograms represented as \code{"histogram"} S3 objects, compute 
#' Wasserstein barycenter. We need one requirement that all histograms in an 
#' input list \code{hists} must have \bold{same breaks}. See the example on how to 
#' construct a histogram on predefined breaks/bins.
#' 
#' @param hists a length-\eqn{N} list of histograms (\code{"histogram"} object) of same breaks.
#' @param p an exponent for the order of the distance (default: 2).
#' @param weights a weight of each image; if \code{NULL} (default), uniform weight is set. Otherwise, 
#' it should be a length-\eqn{N} vector of nonnegative weights. 
#' @param lambda a regularization parameter; if \code{NULL} (default), a paper's suggestion 
#' would be taken, or it should be a nonnegative real number.
#' @param ... extra parameters including \describe{
#' \item{abstol}{stopping criterion for iterations (default: 1e-8).}
#' \item{init.vec}{an initial weight vector (default: uniform weight).}
#' \item{maxiter}{maximum number of iterations (default: 496).}
#' \item{nthread}{number of threads for OpenMP run (default: 1).}
#' \item{print.progress}{a logical to show current iteration (default: \code{TRUE}).}
#' }
#' 
#' @return a \code{"histogram"} object of barycenter.
#' 
#' @examples 
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
#' hh = histbary15B(histxy, maxiter=5)
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' barplot(histxy[[1]]$density, col=rgb(0,0,1,1/4), 
#'         ylim=c(0, 0.75), main="Two Histograms")
#' barplot(histxy[[2]]$density, col=rgb(1,0,0,1/4), 
#'         ylim=c(0, 0.75), add=TRUE)
#' barplot(hh$density, main="Barycenter",
#'         ylim=c(0, 0.75))
#' par(opar)
#' 
#' @seealso \code{\link{bary15B}}
#' 
#' @references 
#' \insertRef{benamou_iterative_2015}{T4transport}
#' 
#' @concept histogram
#' @export
histbary15B <- function(hists, p=2, weights=NULL, lambda=NULL, ...){
  # CHECK THE INPUT
  name.f  = "histbary15B"
  check.f = check_hists(hists, name.f)
  
  # COMPUTE DISTANCE MATRIX
  dxy    = as.matrix(stats::dist(matrix(check.f$midpts,ncol=1))) 
  nhists = length(hists)
  
  # OTHER INFORMATION
  myp = max(1, as.double(p))
  mymarginal = check.f$density
  myweights = valid_multiple_weight(weights, nhists, name.f)
  myweights = myweights/base::sum(myweights)
  if ((length(lambda)==0)&&(is.null(lambda))){
    mylambda = 1/(60/(stats::median(dxy)^myp)) # choice of the paper
  } else {
    mylambda = max(100*.Machine$double.eps, as.double(lambda)) 
  } 
  params = list(...)
  pnames = names(params)
  myiter = max(1, round(ifelse((("maxiter")%in%pnames), params$maxiter, 496)))
  mytol  = max(100*.Machine$double.eps, as.double(ifelse(("abstol"%in%pnames), params$abstol, 1e-8)))
  myshow = as.logical(ifelse(("print.progress"%in%pnames), params$print.progress, FALSE))
  mynthr = max(1, round(ifelse(("nthread"%in%pnames), params$nthread, 1))) # OpenMP Threads
  
  nsupport = base::nrow(dxy)
  if ("init.vec" %in% pnames){
    par_init = as.vector(t(params$init.vec))
    par_init = par_init/base::sum(par_init)
    if ((length(par_init)!=nsupport)||(any(par_init < 0))){
      stop(paste0("* histbary15B : 'init.image' should be of matching size as other images with nonnegative values."))
    }
  } else {
    par_init = rep(1/nsupport, nsupport)
  }
  
  # RUN, WRAP, AND RETURN
  myoutput = routine_bary15B(dxy, mymarginal, myweights, myp, mylambda, myiter, mytol, myshow, par_init, mynthr)
  output   = hists[[1]]
  output$density = as.vector(myoutput)
  output$counts  = round(base::sum(output$counts)*as.vector(myoutput))
  output$xname   = "barycenter"
  return(output)
}
