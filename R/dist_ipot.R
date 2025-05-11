#' Wasserstein Distance via Inexact Proximal Point Method
#' 
#' @description
#' The Inexact Proximal Point Method (IPOT) offers a computationally efficient approach to approximating
#' the Wasserstein distance between two empirical measures by iteratively solving a series of regularized
#' optimal transport problems. This method replaces the entropic regularization used in Sinkhorn's algorithm
#' with a proximal formulation that avoids the explicit use of entropy, thereby mitigating numerical instabilities.
#'
#' Let \eqn{C := \|X_m - Y_n\|^p} be the cost matrix, where \eqn{X_m} and \eqn{Y_n} are the support points of two
#' discrete distributions \eqn{\mu} and \eqn{\nu}, respectively. The IPOT algorithm solves a sequence of optimization problems:
#' \deqn{
#'   \Gamma^{(t+1)} = \arg\min_{\Gamma \in \Pi(\mu, \nu)} \langle \Gamma, C \rangle + \lambda D(\Gamma \| \Gamma^{(t)}),
#' }
#' where \eqn{\lambda > 0} is the proximal regularization parameter and \eqn{D(\cdot \| \cdot)} is the Kullback–Leibler
#' divergence. Each subproblem is solved approximately using a fixed number of inner iterations, making the method inexact.
#'
#' Unlike entropic methods, IPOT does not require \eqn{\lambda \rightarrow 0} for convergence to the unregularized Wasserstein
#' solution. It is therefore more robust to numerical precision issues, especially for small regularization parameters,
#' and provides a closer approximation to the true optimal transport cost with fewer artifacts.
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
#' pmw = paste0("Exact plan\n dist=",round(outw$distance,2))
#' pm1 = paste0("IPOT (lambda=1)\n dist=",round(ipt1$distance,2))
#' pm2 = paste0("IPOT (lambda=10)\n dist=",round(ipt2$distance,2))
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(outw$plan, axes=FALSE, main=pmw)
#' image(ipt1$plan, axes=FALSE, main=pm1)
#' image(ipt2$plan, axes=FALSE, main=pm2)
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{xie_2020_FastProximalPoint}{T4transport}
#' 
#' @concept dist
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

