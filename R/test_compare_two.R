#' test comparison
#' 
#' 
#' 
#' @keywords internal
#' @noRd
compare_two <- function(X, Y, useR = TRUE){
  nx = nrow(X); wx = rep(1/nx, nx)
  ny = nrow(Y); wy = rep(1/ny, ny)
  D2 = util_pairwise_sqdist(X, Y)
  
  # compute the plans
  if (useR){
    return(aux_emd(wx, wy, D2))
  } else {
    return(util_plan_emd_C(wx, wy, D2))
  }
}

# m = 300
# n = 150
# X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
# Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y
# 
# plan_R = compare_two(X, Y, useR=TRUE)
# plan_C = compare_two(X, Y, useR=FALSE)
# 
# par(mfrow=c(1,2))
# image(plan_R, main="lpSolve")
# image(plan_C, main="Bonneel")
# norm(plan_R - plan_C, "F")
# 
# microbenchmark::microbenchmark(
#   plan_R = compare_two(X, Y, useR=TRUE),
#   plan_C = compare_two(X, Y, useR=FALSE),
#   times=5L
# )
