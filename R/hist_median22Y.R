#' Wasserstein Median of Histograms by You et al. (2022)
#' 
#' 
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
  nhists   = length(hists)
  nsupport = base::nrow(dxy)
  
  # weights
  myweights = valid_multiple_weight(weights, nhists, name.f)
  myweights = myweights/base::sum(myweights)
  
  # lambda
  myp = 2.0
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

  # ----------------------------------------------------------------------------
  # COMPUTATION
  # initialize
  hist_old    = par_init
  hist_breaks = as.vector(check.f$midpts)
  
  # iterate
  tmp_dists  = rep(0, nhists)
  tmp_weight = rep(0, nhists)
  for (it in 1:myiter){
    # it-1. update the relative weight
    for (i in 1:nhists){
      # distance / break if too close for one of the histograms
      tmp_histvec  = as.vector(check.f$density[[i]])
      tmp_dists[i] = hist_dist2dis(hist_breaks, hist_old, tmp_histvec)
      if (tmp_dists[i] < 100*.Machine$double.eps){
        output   = hists[[1]]
        output$density = as.vector(tmp_histvec)
        output$counts  = round(base::sum(output$counts)*as.vector(tmp_histvec))
        output$xname   = "barycenter"
        return(output)
      }
      tmp_weight[i] = myweights[i]/as.double(tmp_dists[i])
    }
    tmp_weight = tmp_weight/base::sum(tmp_weight)
    
    # it-2. update a histogram
    hist_new = routine_bary15B(dxy, check.f$density, tmp_weight, 
                               2.0, mylambda, 100, 1e-12, FALSE, 
                               par_init, mynthr)
    
    # it-3. compute the error and update
    increment = hist_dist2dis(hist_breaks, hist_old, hist_new)
    hist_old  = hist_new
    if (increment < mytol){
      break
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

