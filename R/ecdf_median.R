#' Wasserstein Median of Empirical CDFs
#' 
#' Given a collection of empirical cumulative distribution functions 
#' \eqn{F^i (x)} for \eqn{i=1,\ldots,N}, compute the Wasserstein median. This is 
#' obtained by a functional variant of the Weiszfeld algorithm on a set of 
#' quantile functions. 
#' 
#' @param ecdfs a length-\eqn{N} list of \code{"ecdf"} objects by [stats::ecdf()].
#' @param weights a weight of each image; if \code{NULL} (default), uniform weight is set. Otherwise, it should be a length-\eqn{N} vector of nonnegative weights.
#' @param ... extra parameters including \describe{
#' \item{abstol}{stopping criterion for iterations (default: 1e-8).}
#' \item{maxiter}{maximum number of iterations (default: 496).}
#' }
#' 
#' @return an \code{"ecdf"} object of the Wasserstein median.
#' 
#' @examples 
#' \donttest{
#' #----------------------------------------------------------------------
#' #                         Tree Gaussians
#' #
#' # Three Gaussian distributions are parametrized as follows.
#' # Type 1 : (mean, sd) = (-4, 1)
#' # Type 2 : (mean, sd) = ( 0, 1/5)
#' # Type 3 : (mean, sd) = (+6, 1/2)
#' #----------------------------------------------------------------------
#' # GENERATE ECDFs
#' ecdf_list = list()
#' ecdf_list[[1]] = stats::ecdf(stats::rnorm(200, mean=-4, sd=1))
#' ecdf_list[[2]] = stats::ecdf(stats::rnorm(200, mean=+4, sd=0.2))
#' ecdf_list[[3]] = stats::ecdf(stats::rnorm(200, mean=+6, sd=0.5))
#' 
#' # COMPUTE THE MEDIAN
#' emeds = ecdfmed(ecdf_list)
#' 
#' # COMPUTE THE BARYCENTER
#' emean = ecdfbary(ecdf_list)
#' 
#' # QUANTITIES FOR PLOTTING
#' x_grid  = seq(from=-8, to=10, length.out=500)
#' y_type1 = ecdf_list[[1]](x_grid)
#' y_type2 = ecdf_list[[2]](x_grid)
#' y_type3 = ecdf_list[[3]](x_grid)
#' 
#' y_bary = emean(x_grid)
#' y_meds = emeds(x_grid)
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' plot(x_grid, y_bary, lwd=3, col="orange", type="l",
#'      main="Wasserstein Median & Barycenter", 
#'      xlab="x", ylab="Fn(x)", lty=2)
#' lines(x_grid, y_meds, lwd=3, col="blue", lty=2)
#' lines(x_grid, y_type1, col="gray50", lty=3)
#' lines(x_grid, y_type2, col="gray50", lty=3)
#' lines(x_grid, y_type3, col="gray50", lty=3)
#' legend("topleft", legend=c("Median","Barycenter"),
#'         lwd=3, lty=2, col=c("blue","orange"))
#' par(opar)
#' }
#' 
#' @concept ecdf
#' @export
ecdfmed <- function(ecdfs, weights=NULL, ...){
  # --------------------------------------------------------------------------
  # INPUT : EXPLICIT
  # list of ecdfs
  if (!ecdf_check(ecdfs)){
    stop("* ecdfmed : 'ecdfs' should be a list of 'ecdf' objects.")
  }
  
  # weight
  name.f = "ecdfmed"
  weight = valid_weight(weights, length(ecdfs), "weights", name.f)
  
  # extract the data
  dat_ecdf = ecdf_quantiles(ecdfs)
  dat_grid = dat_ecdf$grid
  dat_Finv = dat_ecdf$quantiles
  
  # --------------------------------------------------------------------------
  # INPUT : IMPLICIT
  params = list(...)
  pnames = names(params)
  
  if ("maxiter"%in%pnames){
    par_iter = max(5, round(params$maxiter))
  } else {
    par_iter = 496
  }
  if ("abstol"%in%pnames){
    par_tol = max(100*.Machine$double.eps, as.double(params$abstol))
  } else {
    par_tol = 1e-8
  }
  
  # --------------------------------------------------------------------------
  # COMPUTE
  # parameters
  N = base::nrow(dat_Finv)
  P = base::ncol(dat_Finv)
  
  # initialize
  quantile_old = ecdf_wsum(weight, dat_Finv)
  
  # iterate
  myweight = rep(0, N)
  for (it in 1:par_iter){
    # it-1. update relative weight
    for (n in 1:N){
      dist_n = ecdf_2dist(dat_grid, as.vector(dat_Finv[n,]), quantile_old)
      # stop if matching to one of the data
      if (dist_n < 100*.Machine$double.eps){
        x_med = quantile_old
        y_med = c(0, dat_grid)
        return(stats::stepfun(x_med, y_med))
        break
      }
      myweight[n] = weight[n]/dist_n
    }
    
    # it-2. normalize the weight
    myweight = myweight/base::sum(myweight)
    
    # it-3. update the quantile function
    quantile_new = ecdf_wsum(myweight, dat_Finv)
    
    # it-4. error & update
    increment    = sqrt(sum((quantile_new - quantile_old)^2))
    quantile_old = quantile_new 
    if (increment < par_tol){
      break
    }
  }
  
  # --------------------------------------------------------------------------
  # RETURN
  x_med = quantile_old
  y_med = c(0, dat_grid)
  return(stats::stepfun(x_med, y_med))
}
