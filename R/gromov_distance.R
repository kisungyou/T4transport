#' Gromov-Wasserstein Distance
#' 
#' @description
#' Computes the Gromov-Wasserstein (GW) distance between two metric measure spaces. 
#' Given two distance matrices \eqn{D_X} and \eqn{D_Y} along with their
#' respective marginal distributions, the function solves the GW optimization
#' problem to obtain both the distance value and an associated optimal transport
#' plan.  
#' 
#' The GW distance provides a way to compare datasets that may not lie in the
#' same ambient space by focusing on the intrinsic geometric structure encoded
#' in the pairwise distances. This implementation supports multiple optimization
#' schemes, including majorization–minimization (MM), proximal gradient (PG), and
#' Frank–Wolfe (FW).
#' 
#' @param Dx an \eqn{(M\times M)} distance matrix or a \code{\link[stats]{dist}} object of compatible size.
#' @param Dy an \eqn{(N\times N)} distance matrix or a \code{\link[stats]{dist}} object of compatible size.
#' @param wx a length-\eqn{M} marginal density that sums to \eqn{1}. If \code{NULL} (default), uniform weight is set.
#' @param wy a length-\eqn{N} marginal density that sums to \eqn{1}. If \code{NULL} (default), uniform weight is set.
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of iterations (default: 10).}
#' \item{abstol}{stopping criterion for iterations (default: 1e-6).}
#' \item{method}{optimization method to use; can be one of \code{"mm"}, \code{"pg"}, or \code{"fw"} (default).}
#' }
#' 
#' @return a named list containing\describe{
#' \item{distance}{the computed GW distance value.}
#' \item{plan}{an \eqn{(M\times N)} nonnegative matrix for the optimal transport plan.}
#' }
#' 
#' @examples
#' \dontrun{
#' #-------------------------------------------------------------------
#' #                           Description
#' #
#' # * class 1 : iris dataset (columns 1-4) with perturbations
#' # * class 2 : class 1 rotated randomly in R^4
#' # * class 3 : samples from N((0,0), I)
#' #
#' #  We draw 10 empirical measures from each and compare 
#' #  the regular Wasserstein and GW distance. It is expected that 
#' #  the GW distance between class 1 and class 2 is negligible, 
#' #  while the regular Wasserstein distance is large. For simplicity, 
#' #  limit the cardinalities to 20.
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' set.seed(10)
#' 
#' #  prepare empty lists
#' inputs = vector("list", length=30)
#' 
#' #  generate class 1 and 2
#' iris_mat = as.matrix(iris[sample(1:150,20),1:4])
#' for (i in 1:10){
#'   inputs[[i]] = iris_mat + matrix(rnorm(20*4), ncol=4)
#'   inputs[[i+10]] = inputs[[i]]%*%qr.Q(qr(matrix(runif(16), ncol=4)))
#' }
#' #  generate class 3
#' for (j in 21:30){
#'   inputs[[j]] = matrix(rnorm(20*4), ncol=4)
#' }
#' 
#' ## COMPUTE
#' #  empty arrays
#' dist_RW = array(0, c(30, 30))
#' dist_GW = array(0, c(30, 30))
#' 
#' #  compute pairwise distances
#' for (i in 1:29){
#'   X <- inputs[[i]]
#'   Dx <- stats::dist(X)
#'   for (j in (i+1):30){
#'   Y <- inputs[[j]]
#'   Dy <- stats::dist(Y)
#'   dist_RW[i,j] <- dist_RW[j,i] <- wasserstein(X, Y)$distance
#'   dist_GW[i,j] <- dist_GW[j,i] <- gwdist(Dx, Dy)$distance
#'   }
#' }
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(dist_RW, xaxt="n", yaxt="n", main="Regular Wasserstein distance")
#' image(dist_GW, xaxt="n", yaxt="n", main="Gromov-Wasserstein distance")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{memoli_2011_GromovWassersteinDistances}{T4transport}
#' @concept gromov
#' @export
gwdist <- function(Dx, Dy, wx=NULL, wy=NULL, ...){
  ## INPUTS : EXPLICIT
  if (!valid_gromov_single(Dx)){
    stop("* gwdist: inptu 'Dx' should be a valid distance matrix.")
  }
  if (!valid_gromov_single(Dy)){
    stop("* gwdist: inptu 'Dy' should be a valid distance matrix.")
  }
  input_Dx = as.matrix(Dx)
  input_Dy = as.matrix(Dy)
  input_wx = valid_single_marginal(wx, base::nrow(input_Dx), "gwdist")
  input_wy = valid_single_marginal(wy, base::nrow(input_Dy), "gwdist")
  
  ## INPUTS: IMPLICIT
  params = list(...)
  pnames = names(params)
  
  if ("maxiter" %in% pnames){
    par_iter = max(2, round(params$maxiter))
  } else {
    par_iter = 10
  }
  if ("abstol" %in% pnames){
    par_tol = max(100*.Machine$double.eps, as.double(params$abstol))
  } else {
    par_tol = 1e-6
  }
  if ("method" %in% pnames){
    par_method = params$method
    if (!is.character(par_method)){
      stop("* gwdist: 'method' should be a character string.")
    }
    all_methods = c("pg", "mm", "fw")
    par_method = match.arg(tolower(par_method), all_methods)
  } else {
    par_method = "fw"
  }
  
  ## COMPUTE
  computed = cpp_gwdist(input_Dx, input_Dy, 
                        input_wx, input_wy, 
                        par_iter, par_tol, par_method)
  
  ## RETURN
  output = list()
  output$distance = sqrt(max(0, as.double(computed[["est_cost"]])))
  output$plan = as.matrix(computed[["est_plan"]])
  return(output)
}

# X = as.matrix(iris[sample(1:150,100),1:4])
# Y = X%*%qr.Q(qr(matrix(rnorm(16), ncol=4)))
# Dx = dist(X)
# Dy = dist(Y)
# 
# gwdist(Dx, Dy, method="mm")$distance
# gwdist(Dx, Dy, method="pg")$distance
# gwdist(Dx, Dy, method="fw")$distance
# 
# 
# microbenchmark::microbenchmark(mm_dist = gwdist(Dx, Dy, method="mm")$distance,
#                                pg_dist = gwdist(Dx, Dy, method="pg")$distance,
#                               fw_dist = gwdist(Dx, Dy, method="fw")$distance, 
#                               times=5L)
