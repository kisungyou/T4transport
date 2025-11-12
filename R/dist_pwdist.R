#' Procrustes-Wasserstein Distance
#' 
#' 
#' 
#' 
#' @concept dist
#' @export
pwdist <- function(X, Y, wx=NULL, wy=NULL, ...){
  ## INPUTS: EXPLICIT
  if (is.vector(X)){X = matrix(X, ncol=1)}
  if (is.vector(Y)){Y = matrix(Y, ncol=1)}
  if (!is.matrix(X)){    stop("* pwdist : input 'X' should be a matrix.")  }
  if (!is.matrix(Y)){    stop("* pwdist : input 'Y' should be a matrix.")  }
  if (base::ncol(X)!=base::ncol(Y)){
    stop("* pwdist : input 'X' and 'Y' should be of same dimension.")
  }
  m = base::nrow(X)
  n = base::nrow(Y)
  
  wxname = paste0("'",deparse(substitute(wx)),"'")
  wyname = paste0("'",deparse(substitute(wy)),"'")
  fname  = "pwdist"
  par_wx = valid_single_marginal(wx, m, fname)
  par_wy = valid_single_marginal(wy, n, fname) #valid_weight(wy, n, wyname, fname)
  
  ## INPUTS: IMPLICIT
  params = list(...)
  pnames = names(params)
  
  if ("maxiter" %in% pnames){
    par_iter = max(2, round(params$maxiter))
  } else {
    par_iter = 496
  }
  if ("abstol" %in% pnames){
    par_tol = max(100*.Machine$double.eps, as.double(params$abstol))
  } else {
    par_tol = 1e-10
  }
  
  ## COMPUTE
  compute_cpp   = cpp_pwdist(X, par_wx, Y, par_wy, 
                           par_iter, par_tol)
  
  ## RETURN
  output = list()
  output$distance = sqrt(as.double(compute_cpp[["est_cost"]]))
  output$plan = as.matrix(compute_cpp[["est_plan"]])
  output$alignment = as.matrix(compute_cpp[["est_P"]])
  return(output)
}

# m = sample(100:200, 1)  
# n = sample(100:200, 1)
# X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
# Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y
# Z = Y%*%qr.Q(qr(matrix(rnorm(4), ncol=2)))
# 
# ## COMPUTE WITH DIFFERENT ORDERS
# 
# wasserstein(X, Y)$distance
# wasserstein(X, Z)$distance
# pwdist(X, Y)$distance
# pwdist(X, Z)$distance
# pwdist(Y, Z)$distance