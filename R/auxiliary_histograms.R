## AUXILIARY FUNCTIONS REGARDING THE HISTOGRAMS
#  (01) hist_dist2dis : two discretized one's
#  (02) hist_dist2obj : compute the distance between two histogram's.


# (01) hist_dist2dis ------------------------------------------------------
#' @keywords internal
#' @noRd
hist_dist2dis <- function(vec_x, vec_y1, vec_y2, p=2.0){
  # numerical integration
  n = length(vec_x)
  output  = 0
  vec_ysq = ((vec_y1 - vec_y2)^p)
  for (i in 1:(n-1)){
    output = output + ((vec_ysq[i]+vec_ysq[i+1])/2)*(vec_x[i+1]-vec_x[i])
  }
  return((output^(1/p)))
}

# (02) hist_dist2obj ------------------------------------------------------
#' @keywords internal
#' @noRd
hist_dist2obj <- function(hobj1, hobj2, p=2.0){
  # extract information
  vec_x   = hobj1$mids
  vec_y1  = hobj1$density
  vec_y2  = hobj2$density
  
  # integrate using the (1)
  return(hist_dist2dis(hobj1, hobj2, p=p))
}
