# AUXILIARY FUNCTIONS FOR ECDF
# ecdf_check     : check whether an input is a valid list of ecdfs
# ecdf_quantiles : return the processed data (grid size=1000)
# ecdf_wsum      : compute a weighted sum
# ecdf_2dist     : 2-distance between two discretized quantiles
# ecdf_pdist     : p-distance between two discretized quantiles

# ecdf_check --------------------------------------------------------------
#' @keywords internal
#' @noRd
ecdf_check <- function(ecdfs){
  # should be a list
  if (!is.list(ecdfs)){
    return(FALSE)
  } 
  # should contain ECDFs for all elements in a list
  if (!all(unlist(lapply(ecdfs, inherits, "ecdf"))==TRUE)){
    return(FALSE)
  }
  return(TRUE)
}

# ecdf_quantiles ----------------------------------------------------------
#' @keywords internal
#' @noRd
ecdf_quantiles <- function(ecdfs){
  # parameters
  npts  = 1000
  necdf = length(ecdfs)
  
  # grid
  seq_x = seq(from=0, to=1, length.out=npts)
  
  # quantile functions
  mat_q = array(0,c(necdf, npts))
  for (i in 1:necdf){
    mat_q[i,] = as.vector(stats::quantile(ecdfs[[i]], seq_x))
  }
  
  # return
  return(list(grid=seq_x, quantiles=mat_q))
}


# ecdf_wsum ---------------------------------------------------------------
#' @keywords internal
#' @noRd
ecdf_wsum <- function(weight, ymat){
  n = base::nrow(ymat)
  p = base::ncol(ymat)
  
  rel = weight/base::sum(weight)
  output = rep(0, p)
  for (i in 1:n){
    output = output + rel[i]*as.vector(ymat[i,])
  }
  return(output)
}


# ecdf_2dist --------------------------------------------------------------
#' @keywords internal
#' @noRd
ecdf_2dist <- function(xvec, yvec1, yvec2){
  ydsq = ((yvec1-yvec2)^2)
  n    = length(ydsq)
  return(sqrt(sum(((ydsq[1:(n-1)]+ydsq[2:n])/2)*base::diff(xvec))))
}


# ecdf_pdist --------------------------------------------------------------
#' @keywords internal
#' @noRd
ecdf_pdist <- function(xvec, yvec1, yvec2, p){
  ydsq   = abs(yvec1-yvec2)^p
  n      = length(ydsq)
  output = sum(((ydsq[1:(n-1)]+ydsq[2:n])/2)*base::diff(xvec))
  return(output^(1/p))
}