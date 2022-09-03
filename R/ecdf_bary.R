#' Barycenter of Empirical CDFs
#' 
#' Given a collection of empirical cumulative distribution functions 
#' \eqn{F^i (x)} for \eqn{i=1,\ldots,N}, compute the Wasserstein barycenter 
#' of order 2. This is obtained by taking a weighted average on a set of 
#' corresponding quantile functions. 
#' 
#' @param ecdfs a length-\eqn{N} list of \code{"ecdf"} objects by [stats::ecdf()].
#' @param weights a weight of each image; if \code{NULL} (default), uniform weight is set. Otherwise, it should be a length-\eqn{N} vector of nonnegative weights.
#' @param ... extra parameters including \describe{
#' \item{abstol}{stopping criterion for iterations (default: 1e-8).}
#' \item{maxiter}{maximum number of iterations (default: 496).}
#' }
#' 
#' @return an \code{"ecdf"} object of the Wasserstein barycenter.
#' 
#' @examples
#' #----------------------------------------------------------------------
#' #                         Two Gaussians
#' #
#' # Two Gaussian distributions are parametrized as follows.
#' # Type 1 : (mean, var) = (-4, 1/4)
#' # Type 2 : (mean, var) = (+4, 1/4)
#' #----------------------------------------------------------------------
#' # GENERATE ECDFs
#' ecdf_list = list()
#' ecdf_list[[1]] = stats::ecdf(stats::rnorm(200, mean=-4, sd=0.5))
#' ecdf_list[[2]] = stats::ecdf(stats::rnorm(200, mean=+4, sd=0.5))
#' 
#' # COMPUTE THE BARYCENTER OF EQUAL WEIGHTS
#' emean = ecdfbary(ecdf_list)
#' 
#' # QUANTITIES FOR PLOTTING
#' x_grid  = seq(from=-8, to=8, length.out=100)
#' y_type1 = ecdf_list[[1]](x_grid)
#' y_type2 = ecdf_list[[2]](x_grid)
#' y_bary  = emean(x_grid)
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' plot(x_grid, y_bary, lwd=3, col="red", type="l",
#'      main="Barycenter", xlab="x", ylab="Fn(x)")
#' lines(x_grid, y_type1, col="gray50", lty=3)
#' lines(x_grid, y_type2, col="gray50", lty=3)
#' par(opar)
#' 
#' @concept ecdf
#' @export
ecdfbary <- function(ecdfs, weights=NULL, ...){
  # --------------------------------------------------------------------------
  # INPUT : EXPLICIT
  # list of ecdfs
  if (!ecdf_check(ecdfs)){
    stop("* ecdfbary : 'ecdfs' should be a list of 'ecdf' objects.")
  }
  
  # weight
  name.f = "ecdfbary"
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
  # COMPUTE : WEIGHTED SUM
  fmean = ecdf_wsum(weight, dat_Finv)
  
  # --------------------------------------------------------------------------
  # RETURN
  x_bary = fmean
  y_bary = c(0, dat_grid)
  return(stats::stepfun(x_bary, y_bary))
}