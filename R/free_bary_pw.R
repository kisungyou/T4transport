#' Procrustes-Wasserstein Barycenter
#' 
#' 
#' 
#' @concept free_centroid
#' @export
pwbary <- function(atoms, marginals=NULL, weights=NULL, num_support=100, ...){
  ## INPUT : EXPLICIT
  name.f    = "pwbary"
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
  
  #  compute with the internal routine
  support_out = pwbary_internal(
    X_init = init_measure,
    Xs = par_measures,
    ws = par_marginals,
    rs = par_weights,
    maxiter = par_maxiter,
    abstol = par_abstol
  )
  
  # return
  output = list()
  output[["support"]] = support_out
  output[["weight"]] = rep(1/par_numsupport, par_numsupport)
  return(output)
}




# internal routine for pwbary ---------------------------------------------
# X_init: initial support matrix
# Xs: list of input supports
# ws: list of input weights
# rs: vector of relative weights
# maxiter: maximum number of iterations
# abstol: absolute tolerance
#' @keywords internal
#' @noRd
pwbary_internal <- function(X_init, Xs, ws, rs, maxiter, abstol){
  # initialize
  N = length(Xs)
  X_old = X_init
  X_num = base::nrow(X_old)
  X_weight = rep(1/X_num, X_num)
  
  # iterate
  for (it in seq_len(round(maxiter))){
    # step 1 : compute optimal plan and alignment
    opt_plans = vector("list", length=N)
    opt_align = vector("list", length=N)
    for (n in 1:N){
      compute_now = cpp_pwdist(X_old, X_weight, Xs[[n]], ws[[n]], maxiter, abstol)
      opt_plans[[n]] = as.matrix(compute_now[["est_plan"]])
      opt_align[[n]] = as.matrix(compute_now[["est_P"]])
    }
    
    # step 2 : update the support
    X_new = array(0, c(X_num, ncol(X_old)))
    for (n in 1:N){
      X_new = X_new + rs[n]*(diag(1/X_weight)%*%(opt_plans[[n]]%*%(Xs[[n]]%*%opt_align[[n]])))
    }
    
    # stop
    X_inc = base::norm(X_new-X_old,"F")
    X_old = X_new
    if (X_inc < abstol){
      break
    }
  }
  
  # return
  return(X_old)
}