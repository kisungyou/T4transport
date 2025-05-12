# AUXILIARY FUNCTIONS FOR GENERIC COMPUTATION
#
# aux_emd : solve the POT-style emd problem - emd(a,b,C)





# aux_emd -----------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_emd <- function(a,b,C){
  return(as.matrix(wass_lp(C, 1, a, b)$plan))
}