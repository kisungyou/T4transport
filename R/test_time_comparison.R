#' Test Time Comparison
#' 
#' 
#' @keywords internal
#' @noRd
test_time_comparison <- function(X, Y, lambda=0.01){
  n_x = nrow(X)
  n_y = nrow(Y)
  
  wx = rep(1/n_x, n_x)
  wy = rep(1/n_y, n_y)
  C  = compute_pdist2(X, Y)
  
  start_C = Sys.time()
  output_C = util_plan_emd_C(wx,wy,C)
  time_C = Sys.time() - start_C
  
  start_reg = Sys.time()
  output_reg = util_plan_entropic(wx, wy, C, lambda, maxiter=1000, abstol=1e-12)
  time_reg = Sys.time() - start_reg
  
  output=list()
  output$time = list(emd=time_C, entropic=time_reg)
  output$plan = list(emd=output_C, entropic=output_reg)
  return(output)
}

# par_n = 500
# par_lbd = 0.001
# X = matrix(rnorm(par_n*2), ncol=2)
# Y = matrix(cbind(rnorm(par_n, mean=10, sd=0.01), rnorm(par_n)), ncol=2)
# output = test_time_comparison(X, Y, lambda=par_lbd)
# 
# output$time
# par(mfrow=c(1,2), pty="s")
# image(output$plan$emd, main="Exact LP")
# image(output$plan$entropic, main="Entropic")
