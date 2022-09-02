## AUXILIARY FUNCTIONS FOR GAUSSIAN DISTRIBUTIONS
#
#  gauss_check1d
#  gauss_checknd


# gauss_check1d -----------------------------------------------------------
#' @keywords internal
#' @noRd
gauss_check1d <- function(means, vars){
  cond1 = (length(means)==length(vars))
  cond2 = all(vars>0)
  if (cond1&&cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# gauss_checknd -----------------------------------------------------------
#' @keywords internal
#' @noRd
gauss_checknd <- function(means, vars){
  if (!inherits(means, "matrix")){
    return(FALSE)
  }
  n = base::nrow(means)
  p = base::ncol(vars)
  
  if (length(dim(vars))!=3){
    return(FALSE)
  }
  if (dim(vars)[1]!=p){
    return(FALSE)
  }
  if (dim(vars)[2]!=p){
    return(FALSE)
  }
  if (dim(vars)[3]!=n){
    return(FALSE)
  }
  
  for (i in 1:n){
    tgt = vars[,,i]
    if (!isSymmetric(tgt)){
      return(FALSE)
    }
    if (min(base::eigen(tgt)$value) <= 10*.Machine$double.eps){
      return(FALSE)
    }
  }
  return(TRUE)
}
