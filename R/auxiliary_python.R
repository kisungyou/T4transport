## PYTHON-RELATED FUNCTIONS
#  (1) check_pot



# (1) check_pot -----------------------------------------------------------
#' @keywords internal
#' @noRd
check_pot <- function(fname){
  if (!reticulate::py_module_available("ot")){
    stop(paste0("* ",fname," : Python module is not available. Please use 'install_pot()' first."))
  }
}