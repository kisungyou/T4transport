#' Wasserstein Distance by Inexact Proximal Point Method
#' 
#' Due to high computational cost for linear programming approaches to compute 
#' Wasserstein distance, Cuturi (2013) proposed an entropic regularization 
#' scheme as an efficient approximation to the original problem. This comes with 
#' a regularization parameter \eqn{\lambda > 0} in the term
#' \deqn{\lambda h(\Gamma) = \lambda \sum_{m,n} \Gamma_{m,n} \log (\Gamma_{m,n}).} 
#' IPOT algorithm is known to be relatively robust to the choice of 
#' regularization parameter \eqn{\lambda}. Empirical observation says that 
#' very small number of inner loop iteration like \code{L=1} is sufficient.
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
#' \item{L}{small number of inner loop iterations (default: 1).}
#' }
#' 
#' @return a named list containing\describe{
#' \item{distance}{\eqn{\mathcal{W}_p} distance value}
#' \item{iteration}{the number of iterations it took to converge.}
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
#' n = 30
#' X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
#' Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y
#' 
#' ## COMPARE WITH WASSERSTEIN 
#' outw = wasserstein(X, Y)
#' ipt1 = ipot(X, Y, lambda=1)
#' ipt2 = ipot(X, Y, lambda=10)
#' 
#' ## VISUALIZE : SHOW THE PLAN AND DISTANCE
#' pmw = paste0("wasserstein plan ; dist=",round(outw$distance,2))
#' pm1 = paste0("ipot lbd=1 ; dist=",round(ipt1$distance,2))
#' pm2 = paste0("ipot lbd=10; dist=",round(ipt2$distance,2))
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' image(outw$plan, axes=FALSE, main=pmw)
#' image(ipt1$plan, axes=FALSE, main=pm1)
#' image(ipt2$plan, axes=FALSE, main=pm2)
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{xie_fast_2020}{T4transport}
#' 
#' @concept dist_wass
#' @name ipot
#' @rdname ipot
NULL

#' @rdname ipot
#' @export
ipot <- function(X, Y, p=2, wx=NULL, wy=NULL, lambda=1, ...){
  ## INPUTS : EXPLICIT
  if (is.vector(X)){
    X = matrix(X, ncol=1)
  }
  if (is.vector(Y)){
    Y = matrix(Y, ncol=1)
  }
  if (!is.matrix(X)){    stop("* ipot : input 'X' should be a matrix.")  }
  if (!is.matrix(Y)){    stop("* ipot : input 'Y' should be a matrix.")  }
  if (base::ncol(X)!=base::ncol(Y)){
    stop("* ipot : input 'X' and 'Y' should be of same dimension.")
  }
  m = base::nrow(X)
  n = base::nrow(Y)
  
  wxname = paste0("'",deparse(substitute(wx)),"'")
  wyname = paste0("'",deparse(substitute(wy)),"'")
  fname  = "ipot"
  
  par_wx   = valid_weight(wx, m, wxname, fname)
  par_wy   = valid_weight(wy, n, wyname, fname)
  par_p    = max(1, as.double(p))
  par_D    = as.matrix(compute_pdist2(X, Y))
  par_lbd  = max(sqrt(.Machine$double.eps), as.double(lambda))
  
  ## INPUTS : IMPLICIT
  params = list(...)
  pnames = names(params)
  par_iter  = max(1, round(ifelse((("maxiter")%in%pnames), params$maxiter, 496)))
  par_tol   = max(sqrt(.Machine$double.eps), as.double(ifelse(("abstol"%in%pnames), params$abstol, 1e-10)))
  par_inner = max(1, round(ifelse(("L"%in%pnames), params$L, 1)))
  
  # ## MAIN COMPUTATION
  output = cpp_ipot20(par_wx, par_wy, par_D, par_lbd, par_p, par_iter, par_tol, par_inner)
  return(output)
}

#' @rdname ipot
#' @export
ipotD <- function(D, p=2, wx=NULL, wy=NULL, lambda=1, ...){
  ## INPUTS : EXPLICIT
  name.fun = "ipotD"
  name.D   = paste0("'",deparse(substitute(D)),"'")
  name.wx  = paste0("'",deparse(substitute(wx)),"'")
  name.wy  = paste0("'",deparse(substitute(wy)),"'")
  
  par_D  = valid_distance(D, name.D, name.fun)
  
  m = base::nrow(par_D)
  n = base::ncol(par_D)
  
  par_wx = valid_weight(wx, m, name.wx, name.fun)
  par_wy = valid_weight(wy, n, name.wy, name.fun)
  par_p  = max(1, as.double(p))
  par_lbd  = max(sqrt(.Machine$double.eps), as.double(lambda))
  
  ## INPUTS : IMPLICIT
  params = list(...)
  pnames = names(params)
  par_iter = max(1, round(ifelse((("maxiter")%in%pnames), params$maxiter, 496)))
  par_tol  = max(sqrt(.Machine$double.eps), as.double(ifelse(("abstol"%in%pnames), params$abstol, 1e-10)))
  par_inner = max(1, round(ifelse(("L"%in%pnames), params$L, 1)))
  
  
  # ## MAIN COMPUTATION
  output = cpp_ipot20(par_wx, par_wy, par_D, par_lbd, par_p, par_iter, par_tol, par_inner)
  return(output)
}
