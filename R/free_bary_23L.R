#' Free-Support Barycenter by von Lindheim (2023)
#' 
#' @description
#' For a collection of empirical measures \eqn{\lbrace \mu_k\rbrace_{k=1}^K}, this 
#' function implements the free-support barycenter algorithm introduced by \insertCite{vonlindheim_2023_SimpleApproximativeAlgorithms;textual}{T4transport}.
#' The algorithm takes the first input and its marginal as a reference and performs one-step update of the support. 
#' This version implements `reference` algorithm with \eqn{p=2}.
#' 
#' @param atoms a length-\eqn{K} list where each element is an \eqn{(N_k \times P)} matrix of atoms.
#' @param marginals marginal distributions for empirical measures; if \code{NULL} (default), uniform weights are set for all measures. Otherwise, it should be a length-\eqn{K} list where each element is a length-\eqn{N_i} vector of nonnegative weights that sum to 1.
#' @param weights weights for each individual measure; if \code{NULL} (default), each measure is considered equally. Otherwise, it should be a length-\eqn{K} vector.
#' 
#' @return a list with two elements:
#' \describe{
#'   \item{support}{an \eqn{(N_1 \times P)} matrix of barycenter support points (same number of atoms as the first empirical measure).}
#'   \item{weight}{a length-\eqn{N_1} vector representing barycenter weights (copied from the first marginal).}
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
#' #  The barycenter is computed using the first measure as a reference.
#' #  All measures have uniform weights.
#' #  The barycenter function also considers uniform weights.
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
#' fsbary = rbary23L(myatoms)
#' 
#' ## VISUALIZE
#' #  aligned with CRAN convention
#' opar <- par(no.readonly=TRUE)
#' 
#' #  plot the input measures
#' plot(myatoms[[1]], col="gray90", pch=19, cex=0.5, xlim=c(-6,6), ylim=c(-6,6), 
#'      main="Input Measures", xlab="Dimension 1", ylab="Dimension 2")
#' points(myatoms[[2]], col="gray90", pch=19, cex=0.25)
#' points(myatoms[[3]], col="gray90", pch=19, cex=0.25)
#' points(myatoms[[4]], col="gray90", pch=19, cex=0.25)
#' 
#' #  plot the barycenter
#' points(fsbary$support, col="red", cex=0.5, pch=19)
#' par(opar)
#' }
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @concept free_centroid
#' @export
rbary23L <- function(atoms, marginals=NULL, weights=NULL){
  ## INPUT : EXPLICIT
  name.f    = "rbary23L"
  par_measures  = valid_multiple_measures(atoms, base::ncol(atoms[[1]]), name.f)
  num_atoms     = unlist(lapply(atoms, nrow))
  par_marginals = valid_multiple_marginal(marginals, num_atoms, name.f)
  par_weights   = valid_multiple_weight(weights, length(atoms), name.f)
  
  ## MAJOR COMPUTATION
  #  setup
  N = length(atoms)
  
  #  compute EMD maps
  vec_pi = vector("list", length=N)
  vec_pi[[1]] = diag(par_marginals[[1]])
  for (it in 2:N){
    # compute the squared Euclidean distance
    now_sqdist = as.matrix(compute_pdist2(par_measures[[1]], par_measures[[it]]))
    
    # EMD
    vec_pi[[it]] = aux_emd(par_marginals[[1]], par_marginals[[it]], now_sqdist)
  }
  
  # compute the barycenter
  Y = array(0,c(length(par_marginals[[1]]), ncol(par_measures[[1]])))
  for (it in 1:N){
    Y = Y + par_weights[it]*(diag(1/par_marginals[[1]])%*%(vec_pi[[it]]%*%par_measures[[it]]))
  }
  mu_out = par_marginals[[1]]
  
  # return
  output = list()
  output[["support"]] = Y
  output[["weight"]] = as.vector(mu_out)
  return(output)
}