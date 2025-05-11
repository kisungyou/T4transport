#' Barycenter of Gaussian Distributions in \eqn{\mathbb{R}^p}
#' 
#' Given a collection of \eqn{n}-dimensional Gaussian distributions \eqn{N(\mu_i, \Sigma_i^2)} 
#' for \eqn{i=1,\ldots,n}, compute the Wasserstein barycenter of order 2. 
#' For the barycenter computation of variance components, we use a fixed-point 
#' algorithm by \insertCite{alvarez_2016_FixedpointApproachBarycenters;textual}{T4transport}.
#' 
#' @param means an \eqn{(n\times p)} matrix whose rows are mean vectors.
#' @param vars a \eqn{(p\times p\times n)} array where each slice is covariance matrix.
#' @param weights a weight of each image; if \code{NULL} (default), uniform weight is set. Otherwise, it should be a length-\eqn{n} vector of nonnegative weights.
#' @param ... extra parameters including \describe{
#' \item{abstol}{stopping criterion for iterations (default: 1e-8).}
#' \item{maxiter}{maximum number of iterations (default: 496).}
#' }
#' 
#' @return a named list containing \describe{
#' \item{mean}{a length-\eqn{p} vector for mean of the estimated barycenter distribution.}
#' \item{var}{a \eqn{(p\times p)} matrix for variance of the estimated barycenter distribution.}
#' }
#' 
#' @examples 
#' \donttest{
#' #----------------------------------------------------------------------
#' #                         Two Gaussians in R^2
#' #----------------------------------------------------------------------
#' # GENERATE PARAMETERS
#' # means
#' par_mean = rbind(c(-4,0), c(4,0))
#' 
#' # covariances
#' par_vars = array(0,c(2,2,2))
#' par_vars[,,1] = cbind(c(4,-2),c(-2,4))
#' par_vars[,,2] = cbind(c(4,+2),c(+2,4))
#' 
#' # COMPUTE THE BARYCENTER OF EQUAL WEIGHTS
#' gmean = gaussbarypd(par_mean, par_vars)
#' 
#' # GET COORDINATES FOR DRAWING
#' pt_type1 = gaussvis2d(par_mean[1,], par_vars[,,1])
#' pt_type2 = gaussvis2d(par_mean[2,], par_vars[,,2])
#' pt_gmean = gaussvis2d(gmean$mean, gmean$var)
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' plot(pt_gmean, lwd=2, col="red", type="l",
#'      main="Barycenter", xlab="", ylab="", 
#'      xlim=c(-6,6))
#' lines(pt_type1)
#' lines(pt_type2)
#' par(opar)
#' }
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @seealso [T4transport::gaussbary1d()] for univariate case.
#' @concept gaussian
#' @export
gaussbarypd <- function(means, vars, weights=NULL, ...){
  # --------------------------------------------------------------------------
  # INPUT : EXPLICIT
  # data
  if (!gauss_checknd(means, vars)){
    stop("* gaussbarypd : input 'means' and 'vars' are not valid.")
  }
  N = base::nrow(means)
  P = base::ncol(means)
  
  # weight
  name.f = "gaussbarypd"
  weight = valid_weight(weights, N, "weights", name.f)
  
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
  # compute the mean
  out_mean = rep(0, P)
  for (n in 1:N){
    out_mean = out_mean + as.vector(means[n,])*weight[n]
  }
    
  # compute the variance
  out_var  = as.matrix(gauss_spdbary16A(vars, weight, par_tol, par_iter))
  
  # --------------------------------------------------------------------------
  # RETURN
  return(list(mean=out_mean, var=out_var))
}