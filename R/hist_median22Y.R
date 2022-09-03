#' Wasserstein Median of Histograms by You et al. (2022)
#' 
#' Given multiple histograms represented as \code{"histogram"} S3 objects, compute the 
#' Wasserstein median of order 2. We need one requirement that all histograms in an 
#' input list \code{hists} must have \bold{same breaks}. See the example on how to 
#' construct a histogram on predefined breaks/bins.
#' 
#' @param hists a length-\eqn{N} list of histograms (\code{"histogram"} object) of same breaks.
#' @param weights a weight of each image; if \code{NULL} (default), uniform weight is set. Otherwise, 
#' it should be a length-\eqn{N} vector of nonnegative weights. 
#' @param lambda a regularization parameter; if \code{NULL} (default), a paper's suggestion 
#' would be taken, or it should be a nonnegative real number.
#' @param ... extra parameters including \describe{
#' \item{abstol}{stopping criterion for iterations (default: 1e-8).}
#' \item{init.vec}{an initial weight vector (default: uniform weight).}
#' \item{maxiter}{maximum number of iterations (default: 496).}
#' \item{nthread}{number of threads for OpenMP run (default: 1).}
#' \item{print.progress}{a logical to show current iteration (default: \code{FALSE}).}
#' }
#' 
#' @return a \code{"histogram"} object of the Wasserstein median histogram.
#' 
#' @examples 
#' \donttest{
#' #----------------------------------------------------------------------
#' #                      Binned from Two Gaussians
#' #
#' # EXAMPLE : small example for CRAN for visualization purpose.
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
#' hmean = histbary15B(histxy)
#' hmeds = histmed22Y(histxy)
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' barplot(histxy[[1]]$density, col=rgb(0,0,1,1/4), 
#'         ylim=c(0, 1.05), main="Two Histograms")
#' barplot(histxy[[2]]$density, col=rgb(1,0,0,1/4), 
#'         ylim=c(0, 1.05), add=TRUE)
#' barplot(hmean$density, main="Barycenter",
#'         ylim=c(0, 1.05))
#' barplot(hmeds$density, main="Wasserstein Median",
#'         ylim=c(0, 1.05))
#' par(opar)
#' }
#' 
#' @concept histogram
#' @export
histmed22Y <- function(hists, weights=NULL, lambda=NULL, ...){
  # ----------------------------------------------------------------------------
  # INPUT : EXPLICIT
  # histograms
  name.f  = "histmed22Y"
  check.f = check_hists(hists, name.f)
  
  # grid's pairwise distance
  dxy      = as.matrix(stats::dist(matrix(check.f$midpts,ncol=1))) 
  nimage   = length(hists)
  nsupport = base::nrow(dxy)
  
  # weights
  mypi      = valid_multiple_weight(weights, nimage, name.f)
  mypi      = mypi/base::sum(mypi)
  myweights = mypi
  
  # order & discretized marginal densities
  myp = 2.0
  mymarginal = check.f$density

  # lambda
  if ((length(lambda)==0)&&(is.null(lambda))){
    mylambda = 1/(60/(stats::median(dxy)^myp)) # choice of the paper
  } else {
    mylambda = max(100*.Machine$double.eps, as.double(lambda)) 
  } 
  
  # ----------------------------------------------------------------------------
  # INPUT : IMPLICIT
  params = list(...)
  pnames = names(params)
  
  if ("maxiter"%in%pnames){
    myiter = max(1, round(params$maxiter))
  } else {
    myiter = 496
  }
  if ("abstol"%in%pnames){
    mytol = max(100*.Machine$double.eps, params$abstol)
  } else {
    mytol = 1e-8
  }
  round_digit = ceiling(abs(log10(mytol)))
  if ("nthread"%in%pnames){
    mynthr = max(1, round(params$nthread))
  } else {
    mynthr = 1
  }
  if ("init.vec" %in% pnames){
    par_init = as.vector(params$init.vec)
    par_init = par_init/base::sum(par_init)
    if ((length(par_init)!=nsupport)||(any(par_init < 0))){
      stop(paste0("* histmed22Y : 'init.image' should be of matching size as other images with nonnegative values."))
    }
  } else {
    par_init = rep(1/nsupport, nsupport)
  }
  if ("print.progress"%in%pnames){
    myshow = as.logical(params$print.progress)
  } else {
    myshow = FALSE
  }
  # ----------------------------------------------------------------------------
  # COMPUTATION
  # initialize
  hist_old = par_init
  
  # iterate
  for (it in 1:myiter){
    # it-1. update the relative weights
    for (i in 1:length(myweights)){
      # compute the distance
      sinkhorn_run = cpp_sinkhorn13(as.vector(mymarginal[[i]]), hist_old, dxy, mylambda, myp, 100, mytol)
      
      # break if sufficiently close to one of the data
      if (as.double(sinkhorn_run$distance) < sqrt(.Machine$double.eps)){
        output   = hists[[1]]
        output$density = as.vector(mymarginal[[i]])
        output$counts  = round(base::sum(output$counts)*as.vector(mymarginal[[i]]))
        output$xname   = "Wasserstein median"
        return(output)  
      }
      
      # update weights
      myweights[i] = mypi[i]/as.double(sinkhorn_run$distance) # compute distance by Sinkhorn
    }
    
    # it-2. normalize weights
    myweights = myweights/base::sum(myweights)
    
    # it-3. update histograms
    hist_new = routine_bary15B(dxy, mymarginal, myweights, myp, mylambda, 100, mytol, FALSE, hist_old, mynthr)
    hist_new = as.vector(hist_new)
    
    # it-4. compute the error and update
    increment = max(abs(hist_old-hist_new))
    hist_old  = hist_new
    if (increment < mytol){
      if (myshow){
        print(paste0("* histmed22Y : algorithm terminates at iteration ",it," : increment=",round(increment, round_digit),"."))
      }
      break
    }
    if (myshow){
      print(paste0("* histmed22Y : iteration ",it,"/",myiter," complete : increment=",round(increment, round_digit),"."))
    }
  }
  
  # ----------------------------------------------------------------------------
  # RETURN
  output   = hists[[1]]
  output$density = as.vector(hist_old)
  output$counts  = round(base::sum(output$counts)*as.vector(hist_old))
  output$xname   = "Wasserstein median"
  return(output)
}

