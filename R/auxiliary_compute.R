# AUXILIARY FUNCTIONS FOR GENERIC COMPUTATION
#
# aux_emd : solve the POT-style emd problem - emd(a,b,C)





# aux_emd -----------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_emd <- function(a, b, C){
  # preliminary
  m = length(a); ww_m = matrix(a, ncol=1)
  n = length(b); ww_n = matrix(b, nrow=1)
  ones_m = matrix(rep(1,n), ncol=1)
  ones_n = matrix(rep(1,m), nrow=1)
  
  # CVXR pipeline
  # define a plan
  plan = CVXR::Variable(m,n)
  
  # compute
  wd.obj    <- CVXR::Minimize(CVXR::matrix_trace(t(C)%*%plan))
  wd.const1 <- list(plan >= 0)
  wd.const2 <- list(plan%*%ones_m==ww_m, ones_n%*%plan==ww_n)
  wd.prob   <- CVXR::Problem(wd.obj, c(wd.const1, wd.const2))
  wd.solve  <- CVXR::solve(wd.prob, solver="OSQP")
  
  # return if successful
  if (all(wd.solve$status=="optimal")){ # successful
    gamma <- wd.solve$getValue(plan)
    return(gamma)
  } else { # if not successful, resort to lpSolve
    
    # variable transformation
    c  = as.vector(C)
    A1 = base::kronecker(matrix(1,nrow=1,ncol=n), diag(m))
    A2 = base::kronecker(diag(n), matrix(1,nrow=1,ncol=m))
    A  = rbind(A1, A2)
    
    # problem setup
    f.obj = c
    f.con = A
    f.dir = rep("==",nrow(A))
    f.rhs = c(rep(1/m,m),rep(1/n,n))
    f.sol = (lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs))
    
    # compute and return
    gamma = matrix(f.sol$solution, nrow=m)
    return(gamma)
  }
}