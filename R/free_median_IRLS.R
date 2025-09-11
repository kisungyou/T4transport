#' Free-Support Median by IRLS
#' 
#' @description
#' For a collection of empirical measures \eqn{\lbrace \mu_k\rbrace_{k=1}^K}, 
#' the free-support Wasserstein median, a minimizer to the following 
#' functional
#' \deqn{
#' \mathcal{F}(\nu) = \sum_{k=1}^K w_k \mathcal{W}_2 (\nu, \mu_k ),
#' }
#' is computed using the 
#' generic method of iteratively-reweighted least squares (IRLS) method 
#' according to \insertCite{you_2024_WassersteinMedianProbability;textual}{T4transport}. 
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
#' run_median = rmedIRLS(input_measures, num_support = 50)
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
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @concept free_median
#' @export
rmedIRLS <- function(atoms, marginals=NULL, weights=NULL, num_support=100, ...){
  ## INPUT : EXPLICIT
  name.f    = "rmedIRLS"
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
  
  # temporary flag for printing the updates
  debug_printer = FALSE
  
  ## COMPUTE
  #  initialize
  old_run = list()
  old_measure = aux_ginit(par_measures, par_numsupport)
  old_cost = auxiliary_median_cost(old_measure,
                                   par_measures,
                                   par_marginals,
                                   par_weights)
  
  record_history = rep(0, par_maxiter+1)
  record_history[1] = old_cost

  #  iterate
  for (it in 1:par_maxiter){
    # I need an adjusted weight vector for the median case.
    now_weights = rep(0, length(par_weights))
    for (i in 1:length(par_weights)){
      now_dist = T4transport::wasserstein(X=old_measure, 
                                          Y=par_measures[[i]],
                                          p=2,
                                          wx=rep(1/par_numsupport,par_numsupport),
                                          wy=par_marginals[[i]])$distance
      now_weights[i] = par_weights[i] / max(now_dist, 1e-18)
    }
    now_weights = now_weights/base::sum(now_weights)
    
    # compute a target and the cost
    new_measure = cpp_single_barycenter(
      par_measures,
      par_marginals,
      now_weights,
      old_measure
    )
    new_cost = auxiliary_median_cost(new_measure,
                                   par_measures,
                                   par_marginals,
                                   par_weights)
    record_history[it+1] = new_cost

    # update rule: if it went bad, break
    if (new_cost >= old_cost + 1e-8){
      if (debug_printer){
        print("* break due to increased cost.")
      }
      break
    }
    
    # update rule: if small increment or broken, stop
    inc_update = T4transport::wasserstein(old_measure, new_measure)$distance
    if (!is.finite(inc_update)){
      if (debug_printer){
        print("* break due to broken increment") 
      }
      break
    }
    
    # replace
    rel_cost = abs(old_cost - new_cost)/abs(old_cost)
    old_measure = new_measure
    old_cost = new_cost 
    
    # convergence: stop by the update in support points
    if (inc_update < par_abstol){
      if (debug_printer){
        print(paste0("* rmedIRLS: iteration ",it," - stopped: inc_update"))
      }
      break
    }
    if (rel_cost < par_abstol){
      if (debug_printer){
        print(paste0("* rmedIRLS: iteration ",it," - stopped: rel_cost"))
      }
      break
    }
    
    if (debug_printer){
      print(paste0("* rmedIRLS: iteration ",it," complete: cost=",old_cost))
    }
  }
  
  ## POST-PROCESSING
  #  manipulate the cost history
  vec_history = as.vector(record_history[record_history > 10*.Machine$double.eps])
  
  # return
  output = list()
  output[["support"]] = as.matrix(old_measure)
  output[["weight"]]  = rep(1/par_numsupport, par_numsupport)
  output[["history"]] = vec_history
  return(output)
}



# auxiliary functions for the median --------------------------------------
#' @keywords internal
#' @noRd
auxiliary_median_cost <- function(now_measure,
                                  par_measures,
                                  par_marginals,
                                  par_weights){
  # number of measures
  N = length(par_measures)
  
  # cost 
  cost = 0
  for (n in 1:N){
    # compute the 2-Wasserstein distance
    now_dist = T4transport::wasserstein(X=now_measure,
                                        Y=par_measures[[n]],
                                        p=2,
                                        wy=par_marginals[[n]])
    cost = cost + par_weights[n]*(now_dist$distance)
  }
  return(cost)
}




# # personal example
# # generate "n_good" measures from N(0, I)
# #           "n_bad" measures from N(20, I)
# n_pts = 30
# n_goods = 12
# n_bad = 8
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
# run_median = rmedIRLS(input_measures, num_support = n_pts)
# 
# 
# # plot
# base1 = matrix(rnorm(200*2), ncol=2)
# base2 = matrix(rnorm(50*2), ncol=2); base2[,1] = base2[,1]+20
# base_mat = rbind(base1, base2)
# plot(base_mat, col="grey80", pch=19)
# abline(h=0); abline(v=0)
# points(run_bary$support, col="blue", pch=19)
# points(run_median$support, col="red", pch=19)
