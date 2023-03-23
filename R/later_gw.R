#' Gromov-Wasserstein Distance via Entropic Regularization by 
#' 
#' 
#' 
#' 
#' @keywords internal
#' @noRd
gwdist16P <- function(C1, C2, p1=NULL, p2=NULL, cost=c("KL","quadratic"), ...){
  # --------------------------------------------------------------------------
  # CHECK THE INPUT
  # Explicit
  matC1 = valid_single_similarity(C1, "* gwdist16P : 'C1' is not a valid input.")
  matC2 = valid_single_similarity(C2, "* gwdist16P : 'C2' is not a valid input.")
  
  histp1 = valid_single_marginal(p1, base::nrow(matC1), "gwdist16P")
  histp2 = valid_single_marginal(p2, base::nrow(matC2), "gwdist16P")
  
  # Implicit
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
  if ("lambda.auto"%in%pnames){
    lambda_auto = as.logical(params$lambda.auto)
  } else {
    lambda_auto = TRUE
  }
  if ("lambda.value" %in% pnames){
    lambda_value = max(as.double(params$lambda.value), 100*.Machine$double.eps)
  } else {
    lambda_value = 0.1
  }
  if ("objective" %in% pnames){
    cost_fun = match.arg(tolower(params$objective), c("kl","quad"))
  } else {
    cost_fun = "quad"
  }
  
}

