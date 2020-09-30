#' Wasserstein Distance between Empirical Measures
#' 
#' Given two empirical measures \eqn{\mu, \nu} consisting of \eqn{M} and \eqn{N} observations on \eqn{\mathcal{X}}, \eqn{p}-Wasserstein distance for \eqn{p\geq 1} between two empirical measures 
#' is defined as 
#' \deqn{\mathcal{W}_p (\mu, \nu) = \left( \inf_{\gamma \in \Gamma(\mu, \nu)} \int_{\mathcal{X}\times \mathcal{X}} d(x,y)^p d \gamma(x,y) \right)^{1/p}}
#' where \eqn{\Gamma(\mu, \nu)} denotes the collection of all measures/couplings on \eqn{\mathcal{X}\times \mathcal{X}} 
#' whose marginals are \eqn{\mu} and \eqn{\nu} on the first and second factors, respectively. Please see the section 
#' for detailed description on the usage of the function.
#' 
#' @section Usage \code{wasserstein()} function:
#' We assume empirical measures are defined on the Euclidean space \eqn{\mathcal{X}=\mathbf{R}^d},
#' \deqn{\mu = \sum_{m=1}^M \mu_m \delta_{X_m}\quad\textrm{and}\quad \nu = \sum_{n=1}^N \nu_n \delta_{Y_n}} 
#' and the distance metric used here is standard Euclidean norm \eqn{d(x,y) = \|x-y\|}. Here, the 
#' marginals \eqn{(\mu_1,\mu_2,\ldots,\mu_M)} and \eqn{(\nu_1,\nu_2,\ldots,\nu_N)} correspond to 
#' \code{wx} and \code{wy}, respectively.
#' 
#' @section Using \code{wassersteinD()} function:
#' If other distance measures or underlying spaces are one's interests, we have an option for users to provide 
#' a distance matrix \code{D} rather than vectors, where
#' \deqn{D := D_{M\times N} = d(X_m, Y_n)}
#' for flexible modeling.
#' 
#' @param X an \eqn{(M\times P)} matrix of row observations.
#' @param Y an \eqn{(N\times P)} matrix of row observations.
#' @param D an \eqn{(M\times N)} distance matrix \eqn{d(x_m, y_n)} between two sets of observations.
#' @param p an exponent for the order of the distance (default: 2).
#' @param wx a length-\eqn{M} marginal density that sums to \eqn{1}. If \code{NULL} (default), uniform weight is set.
#' @param wy a length-\eqn{N} marginal density that sums to \eqn{1}. If \code{NULL} (default), uniform weight is set.
#' 
#' @return a named list containing\describe{
#' \item{distance}{\eqn{\mathcal{W}_p} distance value}
#' \item{plan}{an \eqn{(M\times N)} nonnegative matrix for the optimal transport plan.}
#' }
#' @examples 
#' #-------------------------------------------------------------------
#' #  Wasserstein Distance between Samples from Two Bivariate Normal
#' #
#' # * class 1 : samples from Gaussian with mean=(-1, -1)
#' # * class 2 : samples from Gaussian with mean=(+1, +1)
#' #-------------------------------------------------------------------
#' ## SMALL EXAMPLE
#' m = 20
#' n = 10
#' X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
#' Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y
#' 
#' ## COMPUTE WITH DIFFERENT ORDERS
#' out1 = wasserstein(X, Y, p=1)
#' out2 = wasserstein(X, Y, p=2)
#' out5 = wasserstein(X, Y, p=5)
#' 
#' ## VISUALIZE : SHOW THE PLAN AND DISTANCE
#' pm1 = paste0("plan p=1; distance=",round(out1$distance,2))
#' pm2 = paste0("plan p=2; distance=",round(out2$distance,2))
#' pm5 = paste0("plan p=5; distance=",round(out5$distance,2))
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' image(out1$plan, axes=FALSE, main=pm1)
#' image(out2$plan, axes=FALSE, main=pm2)
#' image(out5$plan, axes=FALSE, main=pm5)
#' par(opar)
#' 
#' \dontrun{
#' ## COMPARE WITH ANALYTIC RESULTS
#' #  For two Gaussians with same covariance, their 
#' #  2-Wasserstein distance is known so let's compare !
#' 
#' niter = 5000          # number of iterations
#' vdist = rep(0,niter)
#' for (i in 1:niter){
#'   mm = sample(30:50, 1)
#'   nn = sample(30:50, 1)
#'   
#'   X = matrix(rnorm(mm*2, mean=-1),ncol=2)
#'   Y = matrix(rnorm(nn*2, mean=+1),ncol=2)
#'   
#'   vdist[i] = wasserstein(X, Y, p=2)$distance
#'   if (i%%10 == 0){
#'     print(paste0("iteration ",i,"/", niter," complete.")) 
#'   }
#' }
#' 
#' # Visualize
#' opar <- par(no.readonly=TRUE)
#' hist(vdist, main="Monte Carlo Simulation")
#' abline(v=sqrt(8), lwd=2, col="red")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{peyre_computational_2019}{T4transport}
#' 
#' @concept distance
#' @name wasserstein
#' @rdname wasserstein
NULL


#' @rdname wasserstein
#' @export
wasserstein <- function(X, Y, p=2, wx=NULL, wy=NULL){
  ## CHECK INPUTS
  if (is.vector(X)){
    X = matrix(X, ncol=1)
  }
  if (is.vector(Y)){
    Y = matrix(Y, ncol=1)
  }
  if (!is.matrix(X)){    stop("* wasserstein : input 'X' should be a matrix.")  }
  if (!is.matrix(Y)){    stop("* wasserstein : input 'Y' should be a matrix.")  }
  if (base::ncol(X)!=base::ncol(Y)){
    stop("* wasserstein : input 'X' and 'Y' should be of same dimension.")
  }
  m = base::nrow(X)
  n = base::nrow(Y)
  
  wxname = paste0("'",deparse(substitute(wx)),"'")
  wyname = paste0("'",deparse(substitute(wy)),"'")
  fname  = "wasserstein"
  
  par_wx = valid_weight(wx, m, wxname, fname)
  par_wy = valid_weight(wy, n, wyname, fname)
  par_p  = max(1, as.double(p))
  par_D  = as.matrix(compute_pdist2(X, Y))
  
  output = wass_lp(par_D, par_p, par_wx, par_wy)
  return(output)
}
#' @rdname wasserstein
#' @export
wassersteinD <- function(D, p=2, wx=NULL, wy=NULL){
  ## INPUTS : EXPLICIT
  name.fun = "wassersteinD"
  name.D   = paste0("'",deparse(substitute(D)),"'")
  name.wx  = paste0("'",deparse(substitute(wx)),"'")
  name.wy  = paste0("'",deparse(substitute(wy)),"'")
  
  par_D  = valid_distance(D, name.D, name.fun)
  par_wx = valid_weight(wx, base::nrow(D), name.wx, name.fun)
  par_wy = valid_weight(wy, base::ncol(D), name.wy, name.fun)
  par_p  = max(1, as.double(p))
  
  ## RUN
  output = wass_lp(par_D, par_p, par_wx, par_wy)
  return(output)
}
#' @keywords internal
#' @noRd
wass_lp <- function(dxy, p, wx, wy){
  cxy = (dxy^p)
  m   = nrow(cxy)
  n   = ncol(cxy)
  
  c  = as.vector(cxy)
  A1 = base::kronecker(matrix(1,nrow=1,ncol=n), diag(m))
  A2 = base::kronecker(diag(n), matrix(1,nrow=1,ncol=m))
  A  = rbind(A1, A2)
  
  f.obj = c
  f.con = A
  f.dir = rep("==",nrow(A))
  f.rhs = c(rep(1/m,m),rep(1/n,n))
  f.sol = (lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs))
  
  gamma = matrix(f.sol$solution, nrow=m)
  value = (sum(gamma*cxy)^(1/p))
  
  return(list(distance=value, plan=gamma))
}
