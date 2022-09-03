#' Sampling from a Bivariate Gaussian Distribution for Visualization
#' 
#' This function samples points along the contour of an ellipse represented 
#' by mean and variance parameters for a 2-dimensional Gaussian distribution 
#' to help ease manipulating visualization of the specified distribution. For example, 
#' you can directly use a basic \code{plot()} function directly for drawing.
#' 
#' @param mean a length-\eqn{2} vector for mean parameter.
#' @param var a \eqn{(2\times 2)} matrix for covariance parameter.
#' @param n the number of points to be drawn (default: 500).
#' 
#' @return an \eqn{(n\times 2)} matrix. 
#' 
#' @examples 
#' \donttest{
#' #----------------------------------------------------------------------
#' #                        Three Gaussians in R^2
#' #----------------------------------------------------------------------
#' # MEAN PARAMETERS
#' loc1 = c(-3,0)
#' loc2 = c(0,5)
#' loc3 = c(3,0)
#' 
#' # COVARIANCE PARAMETERS
#' var1 = cbind(c(4,-2),c(-2,4))
#' var2 = diag(c(9,1))
#' var3 = cbind(c(4,2),c(2,4))
#' 
#' # GENERATE POINTS
#' visA = gaussvis2d(loc1, var1)
#' visB = gaussvis2d(loc2, var2)
#' visC = gaussvis2d(loc3, var3)
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' plot(visA[,1], visA[,2], type="l", xlim=c(-5,5), ylim=c(-2,9),
#'      lwd=3, col="red", main="3 Gaussian Distributions")
#' lines(visB[,1], visB[,2], lwd=3, col="blue")
#' lines(visC[,1], visC[,2], lwd=3, col="orange")
#' legend("top", legend=c("Type 1","Type 2","Type 3"),
#'        lwd=3, col=c("red","blue","orange"), horiz=TRUE)
#' par(opar)
#' }
#' 
#' @concept other
#' @export
gaussvis2d <- function(mean, var, n=500){
  # --------------------------------------------------------------------------
  # INPUTS
  # set
  par_mean = as.vector(mean)
  par_vars = as.matrix(var)
  par_npts = max(10, round(n))
  
  # check
  if (length(par_mean)!=2){
    stop("* gaussvis2d : 'mean' should be a vector of length 2.")
  }
  cond1 = (base::nrow(par_vars)==2)
  cond2 = (base::ncol(par_vars)==2)
  cond3 = isSymmetric(par_vars)
  if (!(cond1&&cond2&&cond3)){
    stop("* gaussvis2d : 'var' should be a (2x2) covariance matrix.")
  }
  eig_var = base::eigen(par_vars)
  if (min(eig_var$value) <= .Machine$double.eps){
    stop("* gaussvis2d : 'var' should be of rank 2.")
  }
  
  # extract elements
  a = par_vars[1,1]
  b = par_vars[1,2]
  c = par_vars[2,2]
  
  # --------------------------------------------------------------------------
  # COMPUTE
  # quantity : radii
  lbd1 = ((a+c)/2) + sqrt((((a-c)/2)^2) + (b^2))
  lbd2 = ((a+c)/2) - sqrt((((a-c)/2)^2) + (b^2))
  
  # quantity : rotation angle
  if (b==0){
    if (a>=c){
      theta = 0
    } else {
      theta = pi/2
    }
  } else {
    theta = base::atan2((lbd1-a),b)
  }
  
  # use parametric equation
  vect = seq(from=0, to=2*pi, length.out=par_npts)
  xt = sqrt(lbd1)*cos(theta)*cos(vect) - sqrt(lbd2)*sin(theta)*sin(vect)
  yt = sqrt(lbd1)*sin(theta)*cos(vect) + sqrt(lbd2)*cos(theta)*sin(vect)
  
  # --------------------------------------------------------------------------
  # RETURN
  output = cbind(xt, yt)
  colnames(output) = NULL
  for (i in 1:base::nrow(output)){
    output[i,] = output[i,] + par_mean
  }
  return(output)
}
# draw an ellipse : https://cookierobotics.com/007/