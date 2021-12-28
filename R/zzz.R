.pkgenv <- new.env(parent = emptyenv())

# # RETICULATE : https://cran.r-project.org/web/packages/reticulate/vignettes/package.html
# numpy      <- NULL
# cython     <- NULL
# matplotlib <- NULL
# ot         <- NULL
# 
# .onLoad <- function(libname, pkgname){
#   # use superassignment to update global reference
#   numpy <<- reticulate::import("numpy", delay_load = TRUE)
#   cython <<- reticulate::import("cython", delay_load = TRUE)
#   matplotlib <<- reticulate::import("matplotlib", delay_load = TRUE)
#   ot    <<- reticulate::import("ot",    delay_load = TRUE)
#   reticulate::configure_environment(pkgname)
# }

.onAttach <- function(...){
  ## Retrieve Year Information
  date <- date()
  x <- regexpr("[0-9]{4}", date)
  this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
  
  # Retrieve Current Version
  this.version = packageVersion("T4transport")
  
  ## Print on Screen
  packageStartupMessage("** ------------------------------------------------------- **")
  packageStartupMessage("**   T4transport || Computational Optimal Transport in R ")
  packageStartupMessage("**")
  packageStartupMessage("** Maintainer : Kisung You  (kisungyou@outlook.com)")
  packageStartupMessage("** Version    : ",this.version,"       (",this.year,")",sep="")
  packageStartupMessage("** Website    : https://kisungyou.com/T4transport/")
  packageStartupMessage("**")
  packageStartupMessage("** Please share any bugs or suggestions to the maintainer.")
  packageStartupMessage("** ------------------------------------------------------- **")
}

.onUnload <- function(libpath) {
  library.dynam.unload("T4transport", libpath)
}
