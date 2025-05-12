# AUXILIARY FUNCTIONS FOR GENERIC COMPUTATION
#
# aux_emd : solve the POT-style emd problem - emd(a,b,C)





# aux_emd -----------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_emd <- function(C, a, b){
  return(wass_lp(C, 1, a, b)$plan)
}