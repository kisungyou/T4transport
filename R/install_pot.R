#' Install Python Optimal Transport Library
#' 
#' 
#' @examples 
#' \dontrun{
#' T4transport::install_pot()
#' }
#' 
#' @export
install_pot <- function(method = "auto", conda = "auto"){
  listed = c("numpy","cython","matplotlib","pot")
  for (i in listed){
    if (!reticulate::py_module_available(i)){
      reticulate::py_install(i, method = method, conda = conda)
    }
  }
}