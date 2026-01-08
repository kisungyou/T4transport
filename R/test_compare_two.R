#' test comparison
#' 
#' 
#' 
#' @keywords internal
#' @noRd
compare_two <- function(X, Y, useR = TRUE){
  # nx = nrow(X); wx = rep(1/nx, nx)
  # ny = nrow(Y); wy = rep(1/ny, ny)
  # D2 = util_pairwise_sqdist(X, Y)
  # 
  # # compute the plans
  # output1 = util_plan_emd_C(wx, wy, D2)
  # output2 = util_plan_emd_BSP(wx, wy, D2, nb_plans=20)
  # 
  # result1 = list(plan=output1, distance=sqrt(sum(output1*D2)))
  # result2 = list(plan=output2, distance=sqrt(sum(output2*D2)))
  # results = list(exact=result1, 
  #                approx=result2)
  # return(results)
}
# m = 300
# n = 150
# X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
# Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y
# 
# run_two = compare_two(X, Y)
# 
# par(mfrow=c(1,2))
# image(run_two$exact$plan, main="plan-exact")
# image(run_two$approx$plan, main="plan-approx")
# 
# run_two$approx$distance
# run_two$exact$distance
