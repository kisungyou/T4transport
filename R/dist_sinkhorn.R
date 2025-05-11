#' Wasserstein Distance via Entropic Regularization and Sinkhorn Algorithm
#' 
#' @description
#' To alleviate the computational burden of solving the exact optimal transport problem via linear programming,
#' Cuturi (2013) introduced an entropic regularization scheme that yields a smooth approximation to the
#' Wasserstein distance. Let \eqn{C:=\|X_m - Y_n\|^p} be the cost matrix, where \eqn{X_m} and \eqn{Y_n} are the observations from two distributions \eqn{\mu} and \eqn{nu}. 
#' Then, the regularized problem adds a penalty term to the objective function:
#' \deqn{
#'   W_{p,\lambda}^p(\mu, \nu) = \min_{\Gamma \in \Pi(\mu, \nu)} \langle \Gamma, C \rangle + \lambda \sum_{m,n} \Gamma_{m,n} \log (\Gamma_{m,n}),
#' }
#' where \eqn{\lambda > 0} is the regularization parameter and \eqn{\Gamma} denotes a transport plan.
#' As \eqn{\lambda \rightarrow 0}, the regularized solution converges to the exact Wasserstein solution,
#' but small values of \eqn{\lambda} may cause numerical instability due to underflow.
#' In such cases, the implementation halts with an error; users are advised to increase \eqn{\lambda}
#' to maintain numerical stability.
#' 
#' @param X an \eqn{(M\times P)} matrix of row observations.
#' @param Y an \eqn{(N\times P)} matrix of row observations.
#' @param D an \eqn{(M\times N)} distance matrix \eqn{d(x_m, y_n)} between two sets of observations.
#' @param p an exponent for the order of the distance (default: 2).
#' @param wx a length-\eqn{M} marginal density that sums to \eqn{1}. If \code{NULL} (default), uniform weight is set.
#' @param wy a length-\eqn{N} marginal density that sums to \eqn{1}. If \code{NULL} (default), uniform weight is set.
#' @param lambda a regularization parameter (default: 0.1).
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of iterations (default: 496).}
#' \item{abstol}{stopping criterion for iterations (default: 1e-10).}
#' }
#' 
#' @return a named list containing\describe{
#' \item{distance}{\eqn{\mathcal{W}_p} distance value.}
#' \item{plan}{an \eqn{(M\times N)} nonnegative matrix for the optimal transport plan.}
#' }
#' 
#' @examples 
#' \donttest{
#' #-------------------------------------------------------------------
#' #  Wasserstein Distance between Samples from Two Bivariate Normal
#' #
#' # * class 1 : samples from Gaussian with mean=(-1, -1)
#' # * class 2 : samples from Gaussian with mean=(+1, +1)
#' #-------------------------------------------------------------------
#' ## SMALL EXAMPLE
#' set.seed(100)
#' m = 20
#' n = 10
#' X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
#' Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y
#' 
#' ## COMPARE WITH WASSERSTEIN 
#' outw = wasserstein(X, Y)
#' skh1 = sinkhorn(X, Y, lambda=0.05)
#' skh2 = sinkhorn(X, Y, lambda=0.25)
#' 
#' ## VISUALIZE : SHOW THE PLAN AND DISTANCE
#' pm1 = paste0("Exact Wasserstein:\n distance=",round(outw$distance,2))
#' pm2 = paste0("Sinkhorn (lbd=0.05):\n distance=",round(skh1$distance,2))
#' pm5 = paste0("Sinkhorn (lbd=0.25):\n distance=",round(skh2$distance,2))
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(outw$plan, axes=FALSE, main=pm1)
#' image(skh1$plan, axes=FALSE, main=pm2)
#' image(skh2$plan, axes=FALSE, main=pm5)
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{cuturi_2013_SinkhornDistancesLightspeed}{T4transport}
#' 
#' @concept dist
#' @name sinkhorn
#' @rdname sinkhorn
NULL


#' @rdname sinkhorn
#' @export
sinkhorn <- function(X, Y, p=2, wx=NULL, wy=NULL, lambda=0.1, ...){
  ## INPUTS : EXPLICIT
  if (is.vector(X)){
    X = matrix(X, ncol=1)
  }
  if (is.vector(Y)){
    Y = matrix(Y, ncol=1)
  }
  if (!is.matrix(X)){    stop("* sinkhorn : input 'X' should be a matrix.")  }
  if (!is.matrix(Y)){    stop("* sinkhorn : input 'Y' should be a matrix.")  }
  if (base::ncol(X)!=base::ncol(Y)){
    stop("* sinkhorn : input 'X' and 'Y' should be of same dimension.")
  }
  m = base::nrow(X)
  n = base::nrow(Y)
  
  wxname = paste0("'",deparse(substitute(wx)),"'")
  wyname = paste0("'",deparse(substitute(wy)),"'")
  fname  = "sinkhorn"
  
  par_wx   = valid_weight(wx, m, wxname, fname)
  par_wy   = valid_weight(wy, n, wyname, fname)
  par_p    = max(1, as.double(p))
  par_D    = as.matrix(compute_pdist2(X, Y))
  par_lbd  = max(sqrt(.Machine$double.eps), as.double(lambda))
  
  ## INPUTS : IMPLICIT
  params = list(...)
  pnames = names(params)
  
  if ("maxiter" %in% pnames){
    par_iter = max(2, round(params$maxiter))
  } else {
    par_iter = 496
  }
  if ("abstol" %in% pnames){
    par_tol = max(100*.Machine$double.eps, as.double(params$abstol))
  } else {
    par_tol = 1e-10
  }
  
  # ## MAIN COMPUTATION
  output = cpp_sinkhorn13(par_wx, par_wy, par_D, par_lbd, par_p, par_iter, par_tol)
  return(output)
}


#' @rdname sinkhorn
#' @export
sinkhornD <- function(D, p=2, wx=NULL, wy=NULL, lambda=0.1, ...){
  ## INPUTS : EXPLICIT
  name.fun = "sinkhornD"
  name.D   = paste0("'",deparse(substitute(D)),"'")
  name.wx  = paste0("'",deparse(substitute(wx)),"'")
  name.wy  = paste0("'",deparse(substitute(wy)),"'")
  
  par_D  = valid_distance(D, name.D, name.fun)
  par_wx = valid_weight(wx, base::nrow(D), name.wx, name.fun)
  par_wy = valid_weight(wy, base::ncol(D), name.wy, name.fun)
  par_p  = max(1, as.double(p))
  par_lbd  = max(sqrt(.Machine$double.eps), as.double(lambda))
  
  ## INPUTS : IMPLICIT
  params = list(...)
  pnames = names(params)
  par_iter = max(1, round(ifelse((("maxiter")%in%pnames), params$maxiter, 496)))
  par_tol  = max(sqrt(.Machine$double.eps), as.double(ifelse(("abstol"%in%pnames), params$abstol, 1e-10)))
  
  # ## MAIN COMPUTATION
  output = cpp_sinkhorn13(par_wx, par_wy, par_D, par_lbd, par_p, par_iter, par_tol)
  return(output)
}

# m = sample(100:200, 1)
# n = sample(100:200, 1)
# X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
# Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y
# 
# ## COMPUTE WITH DIFFERENT ORDERS
# 
# sqrt(8)
# wasserstein(X, Y)$distance
# sinkhorn(X, Y, lambda=0.001)$distance
# sinkhorn(X, Y, lambda=0.005)$distance
# sinkhorn(X, Y, lambda=0.01)$distance
# sinkhorn(X, Y, lambda=0.05)$distance
# sinkhorn(X, Y, lambda=0.1)$distance
# sinkhorn(X, Y, lambda=1)$distance
# sinkhorn(X, Y, lambda=10)$distance
# sinkhorn(X, Y, lambda=100)$distance
# sinkhorn(X, Y, lambda=1000)$distance


