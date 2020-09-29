#' Wasserstein Barycenter via Entropic Regularization by Cuturi & Doucet (2014)
#' 
#' Given \eqn{K} empirical measures \eqn{\mu_1, \mu_2, \ldots, \mu_K}, 
#' wasserstein barycenter \eqn{\mu^*} is the solution to the following problem 
#' \deqn{\sum_{k=1}^K \pi_k \mathcal{W}_p^p (\mu, \mu_k)}
#' where \eqn{\pi_k}'s are relative weights of empirical measures. Here we assume 
#' either (1) support atoms in Euclidean space are given, or (2) all pairwise distances between 
#' atoms of the fixed support and empirical measures are given. Here the subgradient 
#' is approximated using the entropic regularization. 
#' 
#' @param support an \eqn{(N\times P)} matrix of rows being atoms for the fixed support.
#' @param measures a length-\eqn{K} list where each element is an \eqn{(N_k \times P)} matrix of atoms.
#' @param distances a length-\eqn{K} list where each element is an \eqn{(N\times N_k)} pairwise distance between atoms of the fixed support and given measures.
#' @param marginals marginal distribution for empirical measures; if \code{NULL} (default), uniform weights are set for all measures. Otherwise, it should be a length-\eqn{K} list where each element is a length-\eqn{N_i} vector of nonnegative weights that sum to 1.
#' @param weights weights for each individual measure; if \code{NULL} (default), each measure is considered equally. Otherwise, it should be a length-\eqn{K} vector.
#' @param lambda regularization parameter (default: 0.1).
#' @param p an exponent for the order of the distance (default: 2).
#' @param ... extra parameters including \describe{
#' \item{abstol}{stopping criterion for iterations (default: 1e-10).}
#' \item{init.vec}{an initial vector (default: uniform weight).}
#' \item{maxiter}{maximum number of iterations (default: 496).}
#' \item{print.progress}{a logical to show current iteration (default: FALSE).}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #     Wasserstein Barycenter for Fixed Atoms with Two Gaussians
#' #
#' # * class 1 : samples from Gaussian with mean=(-4, -4)
#' # * class 2 : samples from Gaussian with mean=(+4, +4)
#' # * target support consists of 7 integer points from -6 to 6,
#' #   where ideally, weight is concentrated near 0 since it's average!
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' #  Empirical Measures
#' set.seed(100)
#' dat1 = matrix(rnorm(sample(20:50, 1)*2, mean=-4, sd=0.5),ncol=2)
#' dat2 = matrix(rnorm(sample(20:50, 1)*2, mean=+4, sd=0.5),ncol=2) 
#' 
#' measures = list()
#' measures[[1]] = dat1
#' measures[[2]] = dat2
#' mydata = rbind(dat1, dat2)
#' 
#' #  Fixed Support
#' support = cbind(seq(from=-8,to=8,by=2),
#'                 seq(from=-8,to=8,by=2))
#' ## COMPUTE
#' comp1 = barysinkhorn14(support, measures, lambda=0.5, maxiter=10)
#' comp2 = barysinkhorn14(support, measures, lambda=1,   maxiter=10)
#' comp3 = barysinkhorn14(support, measures, lambda=5,   maxiter=10)
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' barplot(comp1, main="lambda=0.5")
#' barplot(comp2, main="lambda=1")
#' barplot(comp3, main="lambda=5")
#' par(opar)
#' 
#' @references 
#' \insertRef{cuturi_fast_2014}{T4transport}
#' 
#' @concept barycenter
#' @name barysinkhorn14
#' @rdname barysinkhorn14
NULL

#' @rdname barysinkhorn14
#' @export
barysinkhorn14 <- function(support, measures, marginals=NULL, weights=NULL, lambda=0.1, p=2, ...){
  ## INPUT : EXPLICIT
  name.f    = "barysinkhorn14"
  par_support   = valid_matrixed(support, name.f)
  par_measures  = valid_multiple_measures(measures, ncol(support), name.f)
  num_atoms     = unlist(lapply(measures, nrow))
  par_marginals = valid_multiple_marginal(marginals, num_atoms, name.f)
  par_weights   = valid_multiple_weight(weights, length(measures), name.f)
  par_p    = max(1, as.double(p))
  # input is conventional lbd*h(P); following the notation of the paper (1/lbd)*h(P)
  par_lbd  = 1/max(sqrt(.Machine$double.eps), as.double(lambda)) 
  
  ## INPUT : IMPLICIT
  K = length(measures)
  params = list(...)
  pnames = names(params)
  par_iter = max(1, round(ifelse((("maxiter")%in%pnames), params$maxiter, 496)))
  par_tol  = max(100*.Machine$double.eps, as.double(ifelse(("abstol"%in%pnames), params$abstol, 1e-10)))
  par_show = as.logical(ifelse(("print.progress"%in%pnames), params$print.progress, FALSE))
  nsupport = base::nrow(par_support)
  if ("init.vec" %in% pnames){
    par_init = params$init.vec
    par_init = par_init/base::sum(par_init)
    if ((length(par_init)!=nsupport)||(any(par_init < 0))){
      stop(paste0("* barysinkhorn14 : 'init.vec' should be a vector of length ",nsupport," with nonnegative values."))
    }
  } else {
    par_init = rep(1/nsupport, nsupport)
  }

  ## PREPARE
  par_listdxy = list()
  for (k in 1:K){
    par_listdxy[[k]] = compute_pdist2(par_support, par_measures[[k]])
  }
  
  ## COMPUTE
  output = cpp_barysinkhorn14(par_listdxy, par_marginals, par_weights, par_p, par_lbd, par_iter, par_tol, par_show, par_init)
  return(as.vector(output))
}
#' @rdname barysinkhorn14
#' @export
barysinkhorn14D <- function(distances, marginals=NULL, weights=NULL, lambda=0.1, p=2, ...){
  ## INPUT : EXPLICIT
  name.f      = "barysinkhorn14D"
  par_listdxy = valid_multiple_distance(distances, name.f)
  num_atoms   = unlist(lapply(par_listdxy, ncol))
  par_marginals = valid_multiple_marginal(marginals, num_atoms, name.f)
  par_weights   = valid_multiple_weight(weights, length(par_marginals), name.f)
  par_p    = max(1, as.double(p))
  # input is conventional lbd*h(P); following the notation of the paper (1/lbd)*h(P)
  par_lbd  = 1/max(sqrt(.Machine$double.eps), as.double(lambda)) 
  
  ## INPUT : IMPLICIT
  K = length(par_weights)
  params = list(...)
  pnames = names(params)
  par_iter = max(1, round(ifelse((("maxiter")%in%pnames), params$maxiter, 496)))
  par_tol  = max(100*.Machine$double.eps, as.double(ifelse(("abstol"%in%pnames), params$abstol, 1e-10)))
  par_show = as.logical(ifelse(("print.progress"%in%pnames), params$print.progress, FALSE))
  nsupport = base::nrow(par_listdxy[[1]])
  if ("init.vec" %in% pnames){
    par_init = params$init.vec
    par_init = par_init/base::sum(par_init)
    if ((length(par_init)!=nsupport)||(any(par_init < 0))){
      stop(paste0("* barysinkhorn14D : 'init.vec' should be a vector of length ",nsupport," with nonnegative values."))
    }
  } else {
    par_init = rep(1/nsupport, nsupport)
  }
  
  ## COMPUTE
  output = cpp_barysinkhorn14(par_listdxy, par_marginals, par_weights, par_p, par_lbd, par_iter, par_tol, par_show, par_init)
  return(as.vector(output))
}