#' Tools for Computational Optimal Transport
#' 
#' 
#' @noRd
#' @docType package
#' @name T4transport
#' @aliases T4transport-package
#' @import Rdpack
#' @import CVXR
#' @importFrom lpSolve lp
#' @importFrom CVXR Variable Minimize matrix_trace Problem solve
#' @importFrom stats rnorm median dnorm quantile
#' @importFrom utils packageVersion
#' @importFrom Rcpp evalCpp
#' @useDynLib T4transport
NULL
# pack <- "T4transport"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))