#' Free-Support Barycenter by Riemannian Gradient Descent
#' 
#' @description
#' For a collection of empirical measures \eqn{\lbrace \mu_k\rbrace_{k=1}^K}, 
#' the free-support barycenter of order 2, defined as a minimizer of the following 
#' functional,
#' \deqn{
#' \mathcal{F}(\nu) = \sum_{k=1}^K w_k \mathcal{W}_2^2 (\nu, \mu_k ),
#' }
#' is computed using the Riemannian 
#' gradient descent algorithm. The algorithm is based on the formal Riemannian 
#' geometric view of the 2-Wasserstein space according to \insertCite{otto_2001_GeometryDissipativeEvolution;textual}{T4transport}.
#' 
#' @param atoms a length-\eqn{K} list where each element is an \eqn{(N_k \times P)} matrix of atoms.
#' @param marginals marginal distributions for empirical measures; if \code{NULL} (default), uniform weights are set for all measures. Otherwise, it should be a length-\eqn{K} list where each element is a length-\eqn{N_k} vector of nonnegative weights that sum to 1.
#' @param weights weights for each individual measure; if \code{NULL} (default), each measure is considered equally. Otherwise, it should be a length-\eqn{K} vector.
#' @param num_support the number of support points \eqn{M} for the barycenter (default: 100).
#' @param ... extra parameters including \describe{
#' \item{abstol}{stopping criterion for iterations (default: 1e-6).}
#' \item{maxiter}{maximum number of iterations (default: 10).}
#' }
#' 
#' @return a list with three elements:
#' \describe{
#'   \item{support}{an \eqn{(M \times P)} matrix of barycenter support points.}
#'   \item{weight}{a length-\eqn{M} vector of barycenter weights with all entries being \eqn{1/M}.}
#'   \item{history}{a vector of cost values at each iteration.}
#' }
#' 
#' @examples
#' \donttest{
#' #-------------------------------------------------------------------
#' #     Free-Support Wasserstein Barycenter of Four Gaussians
#' #
#' # * class 1 : samples from Gaussian with mean=(-4, -4)
#' # * class 2 : samples from Gaussian with mean=(+4, +4)
#' # * class 3 : samples from Gaussian with mean=(+4, -4)
#' # * class 4 : samples from Gaussian with mean=(-4, +4)
#' #
#' #  All measures have uniform weights.
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' #  Empirical Measures
#' set.seed(100)
#' unif4 = round(runif(4, 100, 200))
#' dat1 = matrix(rnorm(unif4[1]*2, mean=-4, sd=0.5),ncol=2)
#' dat2 = matrix(rnorm(unif4[2]*2, mean=+4, sd=0.5),ncol=2) 
#' dat3 = cbind(rnorm(unif4[3], mean=+4, sd=0.5), rnorm(unif4[3], mean=-4, sd=0.5))
#' dat4 = cbind(rnorm(unif4[4], mean=-4, sd=0.5), rnorm(unif4[4], mean=+4, sd=0.5))
#' 
#' myatoms = list()
#' myatoms[[1]] = dat1
#' myatoms[[2]] = dat2
#' myatoms[[3]] = dat3
#' myatoms[[4]] = dat4
#' 
#' ## COMPUTE
#' fsbary = rbaryGD(myatoms)
#' 
#' ## VISUALIZE
#' #  aligned with CRAN convention
#' opar <- par(no.readonly=TRUE, mfrow=c(1,2))
#' 
#' #  plot the input measures and the barycenter
#' plot(myatoms[[1]], col="gray90", pch=19, cex=0.5, xlim=c(-6,6), ylim=c(-6,6), 
#'      main="Inputs and Barycenter", xlab="Dimension 1", ylab="Dimension 2")
#' points(myatoms[[2]], col="gray90", pch=19, cex=0.25)
#' points(myatoms[[3]], col="gray90", pch=19, cex=0.25)
#' points(myatoms[[4]], col="gray90", pch=19, cex=0.25)
#' points(fsbary$support, col="red", cex=0.5, pch=19)
#' 
#' #  plot the cost history with only integer ticks
#' plot(seq_along(fsbary$history), fsbary$history, type="b", lwd=2, pch=19,
#'      main="Cost History", xlab="Iteration", ylab="Cost", xaxt='n')
#' axis(1, at=seq_along(fsbary$history))
#' par(opar)
#' }
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @concept free_centroid
#' @export
rbaryGD <- function(atoms, marginals=NULL, weights=NULL, num_support=100, ...){
  ## INPUT : EXPLICIT
  name.f    = "rbaryGD"
  par_measures   = valid_multiple_measures(atoms, base::ncol(atoms[[1]]), name.f)
  num_atoms      = unlist(lapply(atoms, nrow))
  par_marginals  = valid_multiple_marginal(marginals, num_atoms, name.f)
  par_weights    = valid_multiple_weight(weights, length(atoms), name.f)
  par_numsupport = max(2, round(num_support))
  
  ## INPUT : IMPLICIT
  params = list(...)
  pnames = names(params)
  
  if ("maxiter" %in% pnames){
    par_maxiter = max(1, round(params[["maxiter"]]))
  } else {
    par_maxiter = 10
  }
  if ("abstol" %in% pnames){
    par_abstol = max(10*.Machine$double.eps, as.double(params[["abstol"]]))
  } else {
    par_abstol = 1e-6
  }
  
  ## COMPUTE
  #  conditionally
  run_R_init <- TRUE
  
  if (run_R_init){
    # initialize using R
    init_measure = aux_ginit(par_measures, par_numsupport)
    cpprun = cpp_free_bary_gradient_init(
      par_measures,
      par_marginals,
      par_weights,
      par_maxiter,
      par_abstol,
      init_measure
    )
  } else {
    #  run the algorithm
    cpprun = cpp_free_bary_gradient(par_measures, 
                                    par_marginals, 
                                    par_weights,
                                    par_numsupport,
                                    par_maxiter,
                                    par_abstol) 
  }
  
  # manipulate the cost history
  vec_history = as.vector(cpprun$cost_history)
  vec_history = vec_history[vec_history > 1e-18]
  
  # return
  output = list()
  output[["support"]] = as.matrix(cpprun$support)
  output[["weight"]]  = rep(1/par_numsupport, par_numsupport)
  output[["history"]] = vec_history
  return(output)
}