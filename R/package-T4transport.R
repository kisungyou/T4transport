#' Computational Optimal Transport in R
#' 
#' @noRd
#' @name T4transport
#' @aliases T4transport-package
#' @importFrom Rdpack reprompt
#' @import Rdpack
#' @import CVXR
#' @importFrom CVXR Variable Minimize matrix_trace Problem solve
#' @importFrom stats rnorm median dnorm quantile cov
#' @importFrom utils packageVersion
#' @importFrom Rcpp evalCpp
#' @useDynLib T4transport, .registration = TRUE
"_PACKAGE"
# pack <- "T4transport"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))