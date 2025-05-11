#' Barycenter of Gaussian Distributions in \eqn{\mathbb{R}}
#' 
#' Given a collection of Gaussian distributions \eqn{\mathcal{N}(\mu_i, \sigma_i^2)} for \eqn{i=1,\ldots,n}, 
#' compute the Wasserstein barycenter of order 2. For the barycenter computation of 
#' variance components, we use a fixed-point algorithm by \insertCite{alvarez_2016_FixedpointApproachBarycenters;textual}{T4transport}.
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
#' \item{mean}{mean of the estimated barycenter distribution.}
#' \item{var}{variance of the estimated barycenter distribution.}
#' }
#' 
#' @examples
#' #----------------------------------------------------------------------
#' #                         Two Gaussians
#' #
#' # Two Gaussian distributions are parametrized as follows.
#' # Type 1 : (mean, var) = (-4, 1/4)
#' # Type 2 : (mean, var) = (+4, 1/4)
#' #----------------------------------------------------------------------
#' # GENERATE PARAMETERS
#' par_mean = c(-4, 4)
#' par_vars = c(0.25, 0.25)
#' 
#' # COMPUTE THE BARYCENTER OF EQUAL WEIGHTS
#' gmean = gaussbary1d(par_mean, par_vars)
#' 
#' # QUANTITIES FOR PLOTTING
#' x_grid  = seq(from=-6, to=6, length.out=200)
#' y_dist1 = stats::dnorm(x_grid, mean=-4, sd=0.5)
#' y_dist2 = stats::dnorm(x_grid, mean=+4, sd=0.5)
#' y_gmean = stats::dnorm(x_grid, mean=gmean$mean, sd=sqrt(gmean$var)) 
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' plot(x_grid, y_gmean, lwd=2, col="red", type="l",
#'      main="Barycenter", xlab="x", ylab="density")
#' lines(x_grid, y_dist1)
#' lines(x_grid, y_dist2)
#' par(opar)
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @seealso [T4transport::gaussbarypd()] for multivariate case.
#' @concept gaussian
#' @export
gaussbary1d <- function(means, vars, weights=NULL, ...){
  # --------------------------------------------------------------------------
  # INPUT : EXPLICIT
  # data
  dat_mean = as.vector(means)
  dat_vars = as.vector(vars)
  if (!gauss_check1d(dat_mean, dat_vars)){
    stop("* gaussbary1d : input 'means' and 'vars' are not valid.")
  }
  N = length(dat_mean)

  # weight
  name.f = "gaussbary1d"
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

  # compute the mean
  out_mean = 0
  for (i in 1:N){
    out_mean = out_mean + dat_mean[i]*weight[i]
  }

  # compute the variance
  out_var  = as.double(gauss_spdbary16A(array_vars, weight, par_tol, par_iter))
  
  # --------------------------------------------------------------------------
  # RETURN
  return(list(mean=out_mean, var=out_var))
}



# ## personal example : can we do interpolation? : damn yes
# #  parameters
# par_mean = c(-4, 4)
# par_vars = c(0.25, 1)
# 
# # compute
# vec_weight = seq(from=0, to=1, length.out=11)
# vec_interp = list()
# for (i in 1:11){
#   vec_ratio = c(vec_weight[i], 1-vec_weight[i])
#   vec_interp[[i]] = gaussbary1d(par_mean, par_vars, vec_ratio)
# }
# 
# # draw
# par(mfrow=c(1,2))
# xval = seq(from=-6, to=7.5, length.out=1000)
# yvis = c(0,0.81)
# # two distributions
# yval1 = stats::dnorm(xval, mean=par_mean[1], sd=sqrt(par_vars[1]))
# yval2 = stats::dnorm(xval, mean=par_mean[2], sd=sqrt(par_vars[2]))
# plot(xval, yval1, "l", ylim=yvis, main="Two Gaussians")
# lines(xval, yval2)
# 
# # interpolations
# spect10 = RColorBrewer::brewer.pal(11, "Spectral") # color palette
# plot(xval,xval,ylim=yvis,col="white",xlab="x",ylab="density",main="Interpolation")
# for (i in 1:11){
#   gobj = vec_interp[[i]]
#   yval = stats::dnorm(xval, mean=gobj$mean, sd=sqrt(gobj$var))
#   lines(xval, yval, col=spect10[i])
# }