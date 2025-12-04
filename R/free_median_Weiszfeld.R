#' Free-Support Median by Weiszfeld Update with Barycentric Projection
#' 
#' @description
#' For a collection of empirical measures \eqn{\lbrace \mu_k\rbrace_{k=1}^K}, 
#' the free-support Wasserstein median, a minimizer to the following 
#' functional
#' \deqn{
#' \mathcal{F}(\nu) = \sum_{k=1}^K w_k \mathcal{W}_2 (\nu, \mu_k ),
#' }
#' is computed using the OT-adapted version of the Weiszfeld algorithm using 
#' the barycentric projection as a means to recover an optimal displacement map.
#' 
#' @param atoms a length-\eqn{K} list where each element is an \eqn{(N_k \times P)} matrix of atoms.
#' @param marginals marginal distributions for empirical measures; if \code{NULL} (default), uniform weights are set for all measures. Otherwise, it should be a length-\eqn{K} list where each element is a length-\eqn{N_i} vector of nonnegative weights that sum to 1.
#' @param weights weights for each individual measure; if \code{NULL} (default), each measure is considered equally. Otherwise, it should be a length-\eqn{K} vector.
#' @param num_support the number of support points \eqn{M} for the barycenter (default: 100).
#' @param ... extra parameters including \describe{
#' \item{abstol}{stopping criterion for iterations (default: 1e-6).}
#' \item{maxiter}{maximum number of iterations (default: 10).}
#' }
#' 
#' @return a list with three elements:
#' \describe{
#'   \item{support}{an \eqn{(M \times P)} matrix of the Wasserstein median's support points.}
#'   \item{weight}{a length-\eqn{M} vector of median's weights with all entries being \eqn{1/M}.}
#'   \item{history}{a vector of cost values at each iteration.}
#' }
#' 
#' @examples
#' \dontrun{
#' #-------------------------------------------------------------------
#' #     Free-Support Wasserstein Median of Multiple Gaussians
#' #
#' # * class 1 : samples from N((0,0),  Id)
#' # * class 2 : samples from N((20,0), Id)
#' #
#' #  We draw 8 empirical measures of size 50 from class 1, and 
#' #  2 from class 2. All measures have uniform weights.
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' #  8 empirical measures from class 1
#' input_measures = vector("list", length=10L)
#' for (i in 1:8){
#'   input_measures[[i]] = matrix(rnorm(50*2), ncol=2)
#' }
#' for (j in 9:10){
#'   base_draw = matrix(rnorm(50*2), ncol=2)
#'   base_draw[,1] = base_draw[,1] + 20
#'   input_measures[[j]] = base_draw
#' }
#' 
#' ## COMPUTE
#' #  compute the Wasserstein median
#' run_median = rmedWB(input_measures, num_support = 50)
#' #  compute the Wasserstein barycenter
#' run_bary   = rbaryGD(input_measures, num_support = 50)
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' 
#' #  draw the base points of two classes
#' base_1 = matrix(rnorm(80*2), ncol=2)
#' base_2 = matrix(rnorm(20*2), ncol=2)
#' base_2[,1] = base_2[,1] + 20
#' base_mat = rbind(base_1, base_2)
#' plot(base_mat, col="gray80", pch=19)
#' 
#' #  auxiliary information
#' title("estimated barycenter and median")
#' abline(v=0); abline(h=0)
#' 
#' #  draw the barycenter and the median
#' points(run_bary$support, col="red", pch=19)
#' points(run_median$support, col="blue", pch=19)
#' par(opar)
#' }
#' 
#' @concept free_centroid
#' @export
rmedWB <- function(atoms, marginals=NULL, weights=NULL, num_support=100, ...){
  ## INPUT : EXPLICIT
  name.f    = "rmedWB"
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
  #  initialize the measure
  init_measure = aux_ginit(par_measures, par_numsupport)
  cpprun = cpp_free_median_PF(
    par_measures,
    par_marginals,
    par_weights,
    par_maxiter,
    par_abstol,
    init_measure
  )
  
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



# # personal example
# # generate "n_good" measures from N(0, I)
# #           "n_bad" measures from N(20, I)
# n_pts = 10
# n_goods = 7
# n_bad = 3
# input_measures = list()
# for (i in 1:n_goods){
#   input_measures[[i]] = matrix(rnorm(n_pts*2), ncol=2)
# }
# for (i in (n_goods+1):(n_goods+n_bad)){
#   translated = matrix(rnorm(n_pts*2), ncol=2)
#   translated[,1] = translated[,1] + 20
#   input_measures[[i]] = translated
# }
# run_bary = rbaryGD(input_measures, num_support = n_pts)
# run_medIRLS = rmedIRLS(input_measures, num_support = n_pts)
# run_medPF = rmedPF(input_measures, num_support = n_pts)
# 
# 
# # plot
# base1 = matrix(rnorm(200*2), ncol=2)
# base2 = matrix(rnorm(50*2), ncol=2); base2[,1] = base2[,1]+20
# base_mat = rbind(base1, base2)
# plot(base_mat, col="grey80", pch=19)
# abline(h=0); abline(v=0)
# points(run_bary$support, col="blue", pch=19)
# points(run_medIRLS$support, col="red", pch=19)
# points(run_medPF$support, col="green", pch=19)



# # runtime comparison
# # personal example
# # generate "n_good" measures from N(0, I)
# #           "n_bad" measures from N(20, I)
# n_pts = 50
# n_goods = 15
# n_bad = 5
# input_measures = list()
# for (i in 1:n_goods){
#   input_measures[[i]] = matrix(rnorm(n_pts*2), ncol=2)
# }
# for (i in (n_goods+1):(n_goods+n_bad)){
#   translated = matrix(rnorm(n_pts*2), ncol=2)
#   translated[,1] = translated[,1] + 20
#   input_measures[[i]] = translated
# }
# 
# microbenchmark::microbenchmark(
#   run_bary = rbaryGD(input_measures, num_support = n_pts),
#   run_medIRLS = rmedIRLS(input_measures, num_support = n_pts),
#   run_medPF = rmedWB(input_measures, num_support = n_pts),
#   times=3
# )
# 
# 
# run_bary = rbaryGD(input_measures, num_support = n_pts)
# run_medIRLS = rmedIRLS(input_measures, num_support = n_pts)
# run_medPF = rmedPF(input_measures, num_support = n_pts)
# 
# 
# # plot
# base1 = matrix(rnorm(200*2), ncol=2)
# base2 = matrix(rnorm(50*2), ncol=2); base2[,1] = base2[,1]+20
# base_mat = rbind(base1, base2)
# plot(base_mat, col="grey80", pch=19)
# abline(h=0); abline(v=0)
# points(run_bary$support, col="blue", pch=19)
# points(run_medIRLS$support, col="red", pch=19)
# points(run_medPF$support, col="green", pch=19)