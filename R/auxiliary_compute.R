# AUXILIARY FUNCTIONS FOR GENERIC COMPUTATION
#
# aux_emd   : solve the POT-style emd problem - emd(a,b,C)
# aux_ginit : given a list of measures, initialize the barycenter support 





# aux_emd -----------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_emd <- function(a,b,C){
  #return(as.matrix(wass_lp(C, 1, a, b)$plan))
  return(util_plan_emd_C(a,b,C))
}

# aux_ginit ---------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_ginit <- function(list_measures, num_samples){
  # number of measures
  N = length(list_measures)
  p = base::ncol(list_measures[[1]])
  
  # compute the weighted mean
  tmp_mean = rep(0, p)
  tmp_cov  = array(0, c(p,p))
  for (n in 1:N){
    tmp_mean = tmp_mean + as.vector(colMeans(list_measures[[n]]))/N
    tmp_cov  = tmp_cov + stats::cov(list_measures[[n]])/N
  }
  
  # jittering
  tmp_cov = diag(diag(tmp_cov) + (1e-8))
  
  # generate the sample
  return(util_mvrnorm(tmp_mean, tmp_cov, num_samples))
}
