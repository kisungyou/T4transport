#' Procrustes-Wasserstein Distance
#' 
#' @description
#' Given two empirical measures
#' \deqn{\mu = \sum_{m=1}^M \mu_m \delta_{X_m}\quad\textrm{and}\quad \nu = \sum_{n=1}^N \nu_n \delta_{Y_n}} in 
#' \eqn{\mathbb{R}^P}, the Procrustes-Wasserstein (PW) distance is defined as follows:
#' \deqn{
#'   PW_2^2(\mu, \nu) = \min_{Q\in \mathcal{O}(P)} W_2^2(\mu, Q_\# \nu),
#' }
#' where \eqn{\mathcal{O}(P)} is the orthogonal group and \eqn{Q_\#\nu} is the pushforward via \eqn{Q}.
#' 
#' @param X an \eqn{(M\times P)} matrix of row observations.
#' @param Y an \eqn{(N\times P)} matrix of row observations.
#' @param wx a length-\eqn{M} marginal density that sums to \eqn{1}. If \code{NULL} (default), uniform weight is set.
#' @param wy a length-\eqn{N} marginal density that sums to \eqn{1}. If \code{NULL} (default), uniform weight is set.
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of iterations (default: 496).}
#' \item{abstol}{stopping criterion for iterations (default: 1e-10).}
#' }
#' 
#' @return a named list containing\describe{
#' \item{distance}{the computed PW distance value.}
#' \item{plan}{an \eqn{(M\times N)} nonnegative matrix for the optimal transport plan.}
#' \item{alignment}{an optimal alignment matrix of size  \eqn{(P\times P)} in \eqn{\mathcal{O}(P)}.}
#' }
#' 
#' 
#' @examples
#' \dontrun{
#' #-------------------------------------------------------------------
#' #                           Description
#' #
#' # * class 1 : samples from N((0,0),  diag(c(4,1/4)))
#' # * class 2 : samples from N((10,0), diag(c(1/4,4)))
#' # * class 3 : samples from N((10,0), diag(c(1/4,4))) randomly rotated
#' #
#' #  We draw 10 empirical measures from each and compare 
#' #  the regular Wasserstein and PW distance.
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' set.seed(10)
#' 
#' #  prepare empty lists
#' inputs = vector("list", length=30)
#' 
#' #  generate
#' random_rot = qr.Q(qr(matrix(runif(4), ncol=2)))
#' for (i in 1:10){
#'   inputs[[i]] = matrix(rnorm(50*2), ncol=2)
#' }
#' for (j in 11:20){
#'   base_draw = matrix(rnorm(50*2), ncol=2)
#'   base_draw[,1] = base_draw[,1] + 10
#'   
#'   inputs[[j]] = base_draw
#'   inputs[[j+10]] = base_draw%*%random_rot
#' }
#' 
#' ## COMPUTE
#' #  empty arrays
#' dist_RW = array(0, c(30, 30))
#' dist_PW = array(0, c(30, 30))
#' 
#' #  compute pairwise distances
#' for (i in 1:29){
#'   for (j in (i+1):30){
#'   dist_RW[i,j] <- dist_RW[j,i] <- wasserstein(inputs[[i]], inputs[[j]])$distance
#'   dist_PW[i,j] <- dist_PW[j,i] <- pwdist(inputs[[i]], inputs[[j]])$distance
#'   }
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(dist_RW, xaxt="n", yaxt="n", main="Regular Wasserstein distance")
#' image(dist_PW, xaxt="n", yaxt="n", main="PW distance")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{adamo_2025_DepthLookProcrustesWasserstein}{T4transport}
#' 
#' @concept dist
#' @export
pwdist <- function(X, Y, wx=NULL, wy=NULL, ...){
  ## INPUTS: EXPLICIT
  if (is.vector(X)){X = matrix(X, ncol=1)}
  if (is.vector(Y)){Y = matrix(Y, ncol=1)}
  if (!is.matrix(X)){    stop("* pwdist : input 'X' should be a matrix.")  }
  if (!is.matrix(Y)){    stop("* pwdist : input 'Y' should be a matrix.")  }
  if (base::ncol(X)!=base::ncol(Y)){
    stop("* pwdist : input 'X' and 'Y' should be of same dimension.")
  }
  m = base::nrow(X)
  n = base::nrow(Y)
  
  wxname = paste0("'",deparse(substitute(wx)),"'")
  wyname = paste0("'",deparse(substitute(wy)),"'")
  fname  = "pwdist"
  par_wx = valid_single_marginal(wx, m, fname)
  par_wy = valid_single_marginal(wy, n, fname) #valid_weight(wy, n, wyname, fname)
  
  ## INPUTS: IMPLICIT
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
  
  ## COMPUTE
  compute_cpp   = cpp_pwdist(X, par_wx, Y, par_wy, 
                           par_iter, par_tol)
  
  ## RETURN
  output = list()
  output$distance = sqrt(as.double(compute_cpp[["est_cost"]]))
  output$plan = as.matrix(compute_cpp[["est_plan"]])
  output$alignment = as.matrix(compute_cpp[["est_P"]])
  return(output)
}

# m = sample(100:200, 1)  
# n = sample(100:200, 1)
# X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
# Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y
# Z = Y%*%qr.Q(qr(matrix(rnorm(4), ncol=2)))
# 
# ## COMPUTE WITH DIFFERENT ORDERS
# 
# wasserstein(X, Y)$distance
# wasserstein(X, Z)$distance
# pwdist(X, Y)$distance
# pwdist(X, Z)$distance
# pwdist(Y, Z)$distance