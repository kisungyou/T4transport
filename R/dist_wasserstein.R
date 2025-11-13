#' Wasserstein Distance via Linear Programming
#' 
#' @description
#' Given two empirical measures
#' \deqn{\mu = \sum_{m=1}^M \mu_m \delta_{X_m}\quad\textrm{and}\quad \nu = \sum_{n=1}^N \nu_n \delta_{Y_n},}
#' the \eqn{p}-Wasserstein distance for \eqn{p\geq 1} is posited as the following optimization problem
#' \deqn{
#'   W_p^p(\mu, \nu) = \min_{\pi \in \Pi(\mu, \nu)} \sum_{m=1}^M \sum_{n=1}^N \pi_{mn} \|X_m - Y_n\|^p,
#' }
#' where \eqn{\Pi(\mu, \nu)} denotes the set of joint distributions (transport plans) with marginals \eqn{\mu} and \eqn{\nu}. 
#' This function solves the above problem with linear programming, which is a standard approach for 
#' exact computation of the empirical Wasserstein distance. Please see the section 
#' for detailed description on the usage of the function.
#' 
#' @section Using \code{wasserstein()} function:
#' We assume empirical measures are defined on the Euclidean space \eqn{\mathcal{X}=\mathbb{R}^d},
#' \deqn{\mu = \sum_{m=1}^M \mu_m \delta_{X_m}\quad\textrm{and}\quad \nu = \sum_{n=1}^N \nu_n \delta_{Y_n}} 
#' and the distance metric used here is standard Euclidean norm \eqn{d(x,y) = \|x-y\|}. Here, the 
#' marginals \eqn{(\mu_1,\mu_2,\ldots,\mu_M)} and \eqn{(\nu_1,\nu_2,\ldots,\nu_N)} correspond to 
#' \code{wx} and \code{wy}, respectively.
#' 
#' @section Using \code{wassersteinD()} function:
#' If other distance measures or underlying spaces are one's interests, we have an option for users to provide 
#' a distance matrix \code{D} rather than vectors, where
#' \deqn{D := D_{M\times N} = d(X_m, Y_n)}
#' for arbitrary distance metrics beyond the \eqn{\ell_2} norm.
#' 
#' @param X an \eqn{(M\times P)} matrix of row observations.
#' @param Y an \eqn{(N\times P)} matrix of row observations.
#' @param D an \eqn{(M\times N)} distance matrix \eqn{d(x_m, y_n)} between two sets of observations.
#' @param p an exponent for the order of the distance (default: 2).
#' @param wx a length-\eqn{M} marginal density that sums to \eqn{1}. If \code{NULL} (default), uniform weight is set.
#' @param wy a length-\eqn{N} marginal density that sums to \eqn{1}. If \code{NULL} (default), uniform weight is set.
#' 
#' @return a named list containing\describe{
#' \item{distance}{\eqn{\mathcal{W}_p} distance value.}
#' \item{plan}{an \eqn{(M\times N)} nonnegative matrix for the optimal transport plan.}
#' }
#' 
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
#' pm1 = paste0("Order p=1\n distance=",round(out1$distance,2))
#' pm2 = paste0("Order p=2\n distance=",round(out2$distance,2))
#' pm5 = paste0("Order p=5\n distance=",round(out5$distance,2))
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
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
#' niter = 1000          # number of iterations
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
#' \insertRef{peyre_2019_ComputationalOptimalTransport}{T4transport}
#' 
#' @concept dist
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
  
  par_wx = valid_single_marginal(wx, m, fname)
  par_wy = valid_single_marginal(wy, n, fname) #valid_weight(wy, n, wyname, fname)
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
  
  m = base::nrow(par_D)
  n = base::ncol(par_D)
  
  #valid_weight(wy, n, wyname, fname)
  par_wx = valid_single_marginal(wx, m, name.fun)
  par_wy = valid_single_marginal(wy, n, name.fun) 
  par_p  = max(1, as.double(p))
  
  ## RUN
  output = wass_lp(par_D, par_p, par_wx, par_wy)
  return(output)
}
#' @keywords internal
#' @noRd
wass_lp <- function(dxy, p, wx, wy){
  # use the optimized version of Bonneel
  cxy = (dxy^p)
  est_plan = util_plan_emd_C(wx, wy, cxy)
  est_cost = sum(est_plan*cxy)^(1/p)
  return(list(distance=est_cost, plan=est_plan))
}
# wass_lp <- function(dxy, p, wx, wy){
#   # # OLDER VERSION : LPSOLVE
#   # cxy = (dxy^p)
#   # m   = nrow(cxy)
#   # n   = ncol(cxy)
#   # 
#   # c  = as.vector(cxy)
#   # A1 = base::kronecker(matrix(1,nrow=1,ncol=n), diag(m))
#   # A2 = base::kronecker(diag(n), matrix(1,nrow=1,ncol=m))
#   # A  = rbind(A1, A2)
#   # 
#   # f.obj = c
#   # f.con = A
#   # f.dir = rep("==",nrow(A))
#   # f.rhs = c(rep(1/m,m),rep(1/n,n))
#   # f.sol = (lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs))
#   # 
#   # gamma = matrix(f.sol$solution, nrow=m)
#   # value = (sum(gamma*cxy)^(1/p))
#   
#   # NEW VERSION : CVXR
#   # mm = sample(30:50, 1)
#   # nn = sample(30:50, 1)
#   # X = matrix(rnorm(mm*2, mean=-1),ncol=2)
#   # Y = matrix(rnorm(nn*2, mean=+1),ncol=2)
#   # dxy = array(0,c(mm,nn))
#   # for (i in 1:mm){
#   #   for (j in 1:nn){
#   #     dxy[i,j] <- sqrt(sum((as.vector(X[i,])-as.vector(Y[j,]))^2))
#   #   }
#   # }
#   # wx = rep(1/mm, mm)
#   # wy = rep(1/nn, nn)
#   # p  = 2
#   
#   cxy = (dxy^p)
#   m   = length(wx); ww_m = matrix(wx, ncol=1)
#   n   = length(wy); ww_n = matrix(wy, nrow=1)
#   ones_m = matrix(rep(1,n),ncol=1)
#   ones_n = matrix(rep(1,m),nrow=1)
#   plan   = CVXR::Variable(m,n)
# 
#   #wd.obj    <- CVXR::Minimize(CVXR::matrix_trace(t(cxy)%*%plan))
#   wd.obj    <- CVXR::Minimize(CVXR::sum_entries(CVXR::multiply(cxy, plan)))
#   wd.const1 <- list(plan >= 0)
#   wd.const2 <- list(plan%*%ones_m==ww_m, ones_n%*%plan==ww_n)
#   wd.prob   <- CVXR::Problem(wd.obj, c(wd.const1, wd.const2))
#   wd.solve  <- CVXR::solve(wd.prob, solver="OSQP")
#   
#   if (all(wd.solve$status=="optimal")){ # successful
#     gamma <- wd.solve$getValue(plan)
#     value <- (base::sum(gamma*cxy)^(1/p))
#     
#     return(list(distance=value, plan=gamma))
#   } else {                              # failed : use lpsolve
#     cxy = (dxy^p)
#     m   = nrow(cxy)
#     n   = ncol(cxy)
#     
#     c  = as.vector(cxy)
#     A1 = base::kronecker(matrix(1,nrow=1,ncol=n), diag(m))
#     A2 = base::kronecker(diag(n), matrix(1,nrow=1,ncol=m))
#     A  = rbind(A1, A2)
#     
#     f.obj = c
#     f.con = A
#     f.dir = rep("==",nrow(A))
#     f.rhs = c(rep(1/m,m),rep(1/n,n))
#     f.sol = (lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs))
#     
#     gamma = matrix(f.sol$solution, nrow=m)
#     value = (sum(gamma*cxy)^(1/p))
#     
#     return(list(distance=value, plan=gamma))
#   }
# }
