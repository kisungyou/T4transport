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
  output1 = util_plan_emd_C(wx, wy, D2)
  output2 = util_dual_emd_C(wx, wy, D2, TRUE)
  return(list(plan=output1,
              dual=output2))
}
# m = 300
# n = 150
# X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
# Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y
# 
# run_two = compare_two(X, Y)
# 
# par(mfrow=c(1,2))
# image(run_two$plan, main="plan")
# image(run_two$dual$G, main="dual plan")
