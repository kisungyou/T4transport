#' Wasserstein Median of Gaussian Distributions in \eqn{\mathbf{R}}
#' 
#' Given a collection of Gaussian distributions \eqn{\mathcal{N}(\mu_i, \sigma_i^2)} for \eqn{i=1,\ldots,n}, 
#' compute the Wasserstein median.
#' 
#' @param means a length-\eqn{n} vector of mean parameters.
#' @param vars a length-\eqn{n} vector of variance parameters.
#' @param weights a weight of each image; if \code{NULL} (default), uniform weight is set. Otherwise, it should be a length-\eqn{n} vector of nonnegative weights.
#' @param ... extra parameters including \describe{
#' \item{abstol}{stopping criterion for iterations (default: 1e-8).}
#' \item{maxiter}{maximum number of iterations (default: 496).}
#' }
#' 
#' @return a named list containing \describe{
#' \item{mean}{mean of the estimated median distribution.}
#' \item{var}{variance of the estimated median distribution.}
#' }
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
#' # GENERATE PARAMETERS
#' par_mean = c(-4, 0, +6)
#' par_vars = c(1, 0.04, 0.25)
#' 
#' # COMPUTE THE WASSERSTEIN MEDIAN
#' gmeds = gaussmed1d(par_mean, par_vars)
#' 
#' # COMPUTE THE BARYCENTER 
#' gmean = gaussbary1d(par_mean, par_vars)
#' 
#' # QUANTITIES FOR PLOTTING
#' x_grid  = seq(from=-6, to=8, length.out=1000)
#' y_dist1 = stats::dnorm(x_grid, mean=par_mean[1], sd=sqrt(par_vars[1]))
#' y_dist2 = stats::dnorm(x_grid, mean=par_mean[2], sd=sqrt(par_vars[2]))
#' y_dist3 = stats::dnorm(x_grid, mean=par_mean[3], sd=sqrt(par_vars[3]))
#' 
#' y_gmean = stats::dnorm(x_grid, mean=gmean$mean, sd=sqrt(gmean$var)) 
#' y_gmeds = stats::dnorm(x_grid, mean=gmeds$mean, sd=sqrt(gmeds$var))
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' plot(x_grid, y_gmeds, lwd=3, col="red", type="l",
#'      main="Three Gaussians", xlab="x", ylab="density", 
#'      xlim=range(x_grid), ylim=c(0,2.5))
#' lines(x_grid, y_gmean, lwd=3, col="blue")
#' lines(x_grid, y_dist1, lwd=1.5, lty=2)
#' lines(x_grid, y_dist2, lwd=1.5, lty=2)
#' lines(x_grid, y_dist3, lwd=1.5, lty=2)
#' legend("topleft", legend=c("Median","Barycenter"),
#'        col=c("red","blue"), lwd=c(3,3), lty=c(1,2))
#' par(opar)
#' }
#' 
#' @seealso [T4transport::gaussmedpd()] for multivariate case.
#' @concept gaussian
#' @export
gaussmed1d <- function(means, vars, weights=NULL, ...){
  # --------------------------------------------------------------------------
  # INPUT : EXPLICIT
  # data
  dat_mean = as.vector(means)
  dat_vars = as.vector(vars)
  if (!gauss_check1d(dat_mean, dat_vars)){
    stop("* gaussmed1d : input 'means' and 'vars' are not valid.")
  }
  N = length(dat_mean)
  
  # weight
  name.f = "gaussmed1d"
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
  # re-arrange the data
  array_mean = matrix(dat_mean, ncol=1)
  array_vars = array(0,c(1,1,N))
  for (i in 1:N){
    array_vars[,,i] = dat_vars[i]
  }
  
  # compute the median
  out_mean = as.double(gauss_weiszfeld(array_mean, weight, par_tol, par_iter));
  
  # compute the variance
  out_var = as.double(gauss_spdmed22Y(array_vars, weight, par_tol, par_iter))

  # --------------------------------------------------------------------------
  # RETURN
  return(list(mean=out_mean, var=out_var))
}

