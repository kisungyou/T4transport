#' Gromov-Wasserstein Barycenter
#' 
#' @description
#' Computes the Gromov–Wasserstein (GW) barycenter of a collection of metric
#' measure spaces. Given a list of distance matrices \eqn{{D^{(k)}}{k=1}^K}
#' and their corresponding marginal distributions, the function estimates a
#' synthetic metric space whose intrinsic geometry best represents the input
#' collection under the GW criterion.
#'
#' The GW barycenter is defined as the minimizer of a multi-measure
#' Gromov–Wasserstein objective, where each dataset contributes according to a
#' user-specified barycentric weight. Since the problem is jointly non-convex in
#' the barycenter metric and the coupling matrices, the algorithm proceeds
#' through an outer–inner iterative procedure.
#' 
#' @param distances a length-\eqn{K} list where each element is either an \eqn{(N_k \times N_k)} distance matrix or an object of class \code{dist} representing the pairwise distances for each empirical measure.
#' @param marginals marginal distributions for empirical measures; if \code{NULL} (default), uniform weights are set for all measures. Otherwise, it should be a length-\eqn{K} list where each element is a length-\eqn{N_k} vector of nonnegative weights that sum to 1.
#' @param weights weights for each individual measure; if \code{NULL} (default), each measure is considered equally. Otherwise, it should be a length-\eqn{K} vector.
#' @param num_support the number of support points \eqn{M} for the barycenter (default: 100).
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of iterations (default: 10).}
#' \item{abstol}{stopping criterion for iterations (default: 1e-6).}
#' \item{method}{optimization method to use; can be one of \code{"mm"}, \code{"pg"}, or \code{"fw"} (default).}
#' }
#' 
#' @return A named list containing \describe{
#' \item{dist}{an object of class \code{dist} representing the GW barycenter.}
#' \item{weight}{a length-\eqn{M} vector of barycenter weights with all entries being \eqn{1/M}.}
#' }
#' 
#' @examples
#' \dontrun{
#' #-------------------------------------------------------------------
#' #                           Description
#' #
#' # GW barycenter computation is quite expensive. In this example, 
#' # we draw a small set of empirical measures from the digit '3'
#' # images and compute their GW barycenter with a small number of 
#' # support points. The attained barycenter distance matrix is then 
#' # passed onto the classical MDS algorithm for visualization.
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' data(digits)
#' data_D = vector("list", length=5)
#' data_W = vector("list", length=5)
#' for (i in 1:5){
#'   img_now = img2measure(digits3[[i]])
#'   data_D[[i]] = stats::dist(img_now$support)
#'   data_W[[i]] = as.vector(img_now$weight)
#' }
#' 
#' ## COMPUTE
#' bary_dist <- gwbary(data_D, marginals=data_W, num_support=100)
#' bary_cmd2 <- stats::cmdscale(bary_dist$dist, k=2)
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(pty="s")
#' plot(bary_cmd2, main="GW Barycenter Embedding",
#'      xaxt="n", yaxt="n", pch=19, xlab="", ylab="")
#' par(opar)
#' }
#' 
#' @concept gromov
#' @export
gwbary <- function(distances, marginals=NULL, weights=NULL, num_support=100, ...){
  ## INPUT : EXPLICIT
  name.f    = "gwbary"
  if (!valid_gromov_list(distances)){
    stop("* gwbary: input 'distances' should be a list of valid distance matrices.")
  }
  par_distances = vector("list", length=length(distances))
  for (i in 1:length(distances)){
    par_distances[[i]] = as.matrix(distances[[i]])
  }
  num_atoms      = unlist(lapply(par_distances, nrow))
  par_marginals  = valid_multiple_marginal(marginals, num_atoms, name.f)
  par_weights    = valid_multiple_weight(weights, length(num_atoms), name.f)
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
  if ("method" %in% pnames){
    par_method = params$method
    if (!is.character(par_method)){
      stop("* gwbary: 'method' should be a character string.")
    }
    all_methods = c("pg", "mm", "fw")
    par_method = match.arg(tolower(par_method), all_methods)
  } else {
    par_method = "fw"
  }
  
  ## COMPUTE
  init_D   <- gwbary_init(par_distances, par_numsupport) # initialization
  computed <- cpp_gwbary(par_distances, 
                         par_marginals, 
                         par_weights,
                         init_D,
                         par_maxiter,
                         par_abstol,
                         par_method)

  # RETURN
  output = list()
  output[["dist"]] = stats::as.dist(computed)
  output[["weight"]]  = rep(1/par_numsupport, par_numsupport)
  return(output)
}

#' @keywords internal
#' @noRd
gwbary_init <- function(list_D, num_support){
  K = length(list_D)
  
  # run cmds
  list_cmds <- vector("list", length=K)
  for (k in 1:K){
    list_cmds[[k]] <- util_cmds(list_D[[k]], 2)
  }
  if (any(unlist(lapply(list_cmds, ncol)) > 2)){
    for (k in 1:K){
      if (ncol(list_cmds[[k]]) < 2){
        vec1 = as.vector(list_cmds[[k]])
        vec2 = stats::rnorm(length(vec1))
        list_cmds[[k]] = cbind(vec1, vec2)
      }
    }
  }
  
  # ginit trick
  sam_measure <- aux_ginit(list_cmds, num_support)
  sam_dist <- as.matrix(stats::dist(sam_measure))
  return(sam_dist)
}