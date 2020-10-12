#' Gromov-Wasserstein Distance
#' 
#' 
#' 
#' 
#' 
#' @return a named list containing\describe{
#' \item{distance}{Gromov-Wasserstein distance value.}
#' \item{plan}{an \eqn{(M\times N)} nonnegative matrix for the optimal transport plan.}
#' }
#' 
#' @examples 
#' #-------------------------------------------------------------------
#' #  Gromow-Wasserstein Distance between Two Displaced Bivariate Normal
#' #
#' # * class 1 : samples from Gaussian with mean=(-5, -5)
#' # * class 2 : samples from Gaussian with mean=(+5, +5)
#' #-------------------------------------------------------------------
#' if (reticulate::py_module_available("ot")){
#'   ## SMALL EXAMPLE
#'   m = sample(50:100, 1)
#'   n = sample(50:100, 1)
#'   X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
#'   Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y
#' 
#'   xyloc = rbind(X, Y)
#'   xycol = c(rep(2,m),rep(4,n))
#' 
#'   ## COMPUTE WITH DIFFERENT COST
#'   out.sq = gw(X, Y, loss="Square")
#'   out.kl = gw(X, Y, loss="kl")
#' 
#'   ## VISUALIZE : SHOW THE PLAN AND DISTANCE
#'   pm.sq = paste0("Square loss; dist=",round(out.sq$distance,4))
#'   pm.kl = paste0("KL loss; dist=",round(out.kl$distance,4))
#' 
#'   opar <- par(no.readonly=TRUE)
#'   par(mfrow=c(1,3))
#'   plot(xyloc, col=xycol, pch=19, main="data")
#'   image(out.sq$plan, axes=FALSE, main=pm.sq)
#'   image(out.kl$plan, axes=FALSE, main=pm.kl)
#'   par(opar)
#' } else {
#'   print("No POT module detected. Use 'install_pot()' to run this example.")
#' }
#' 
#' @references 
#' \insertRef{memoli_gromov_2011}{T4transport}
#' 
#' @concept gromov
#' @rdname gw
#' @name gw
NULL

#' @rdname gw
#' @export
gw <- function(X, Y, wx=NULL, wy=NULL, loss=c("square","KL"), armijo=FALSE){
  ## CHECK OT MODULE AVAILABILITY 
  fname = "gw"
  # check_pot(fname)
  
  ## INPUTS : EXPLICIT
  if (is.vector(X)){
    X = matrix(X, ncol=1)
  }
  if (is.vector(Y)){
    Y = matrix(Y, ncol=1)
  }
  if (!is.matrix(X)){    stop("* gw : input 'X' should be a matrix.")  }
  if (!is.matrix(Y)){    stop("* gw : input 'Y' should be a matrix.")  }
  
  D1 = as.matrix(stats::dist(X)); M = base::nrow(D1)
  D2 = as.matrix(stats::dist(Y)); N = base::nrow(D2)
  new_p1 = matrix(valid_single_marginal(wx, M, fname), ncol=1)
  new_p2 = matrix(valid_single_marginal(wy, N, fname), ncol=1)
  
  # # INPUT : INTERNAL
  # params = list(...)
  # pnames = names(params)
  # myiter = as.integer(max(1, round(ifelse((("maxiter")%in%pnames), params$maxiter, 100))))
  # mytol  = as.double(max(100*.Machine$double.eps, as.double(ifelse(("abstol"%in%pnames), params$abstol, 1e-8))))
  myloss   = match.arg(tolower(loss), c("square","kl"))
  myarmijo = as.logical(armijo)
  
  # PYTHON : WRAP THE DATA
  pCx = reticulate::r_to_py(D1)
  pCy = reticulate::r_to_py(D2)
  ppx = reticulate::r_to_py(new_p1); ppx = ppx$flatten()
  ppy = reticulate::r_to_py(new_p2); ppy = ppy$flatten()
  
  # PYTHON : RUN
  if (all(myloss=="square")){
    pyrun = ot$gromov$gromov_wasserstein(pCx,pCy,ppx,ppy,"square_loss", armijo=myarmijo, verbose=FALSE, log=TRUE)
  } else {
    pyrun = ot$gromov$gromov_wasserstein(pCx,pCy,ppx,ppy,"kl_loss", armijo=myarmijo, verbose=FALSE, log=TRUE)
  }
  
  # RETURN OUTPUT 
  output = list()
  output$plan = pyrun[[1]]
  output$distance = pyrun[[2]]$gw_dist
  return(output)
}

#' @rdname gw
#' @export
gwD <- function(distX, distY, wx=NULL, wy=NULL, loss=c("square","KL"), armijo=FALSE){
  # CHECK OT MODULE AVAILABILITY 
  fname = "gwD"
  # check_pot(fname)
  
  # INPUT : EXPLICIT
  if (!inherits(distX,"dist")){    stop(paste0("* ",fname," : input 'distX' should be a 'dist' object. See help with 'help(dist)'."))  }
  if (!inherits(distY,"dist")){    stop(paste0("* ",fname," : input 'distY' should be a 'dist' object. See help with 'help(dist)'."))  }
  D1  = as.matrix(distX); M = base::nrow(D1)
  D2  = as.matrix(distY); N = base::nrow(D2)

  new_p1 = matrix(valid_single_marginal(wx, M, fname), ncol=1)
  new_p2 = matrix(valid_single_marginal(wy, N, fname), ncol=1)
  
  # # INPUT : INTERNAL
  # params = list(...)
  # pnames = names(params)
  # myiter = as.integer(max(1, round(ifelse((("maxiter")%in%pnames), params$maxiter, 100))))
  # mytol  = as.double(max(100*.Machine$double.eps, as.double(ifelse(("abstol"%in%pnames), params$abstol, 1e-8))))
  myloss   = match.arg(tolower(loss), c("square","kl"))
  myarmijo = as.logical(armijo)

  # PYTHON : WRAP THE DATA
  pCx = reticulate::r_to_py(D1)
  pCy = reticulate::r_to_py(D2)
  ppx = reticulate::r_to_py(new_p1); ppx = ppx$flatten()
  ppy = reticulate::r_to_py(new_p2); ppy = ppy$flatten()
  
  # PYTHON : RUN
  if (all(myloss=="square")){
    pyrun = ot$gromov$gromov_wasserstein(pCx,pCy,ppx,ppy,"square_loss", armijo=myarmijo, verbose=FALSE, log=TRUE)
  } else {
    pyrun = ot$gromov$gromov_wasserstein(pCx,pCy,ppx,ppy,"kl_loss", armijo=myarmijo, verbose=FALSE, log=TRUE)
  }
  
  # RETURN OUTPUT 
  output = list()
  output$plan = pyrun[[1]]
  output$distance = pyrun[[2]]$gw_dist
  return(output)
}

# m = 40
# n = 30
# X = dist(matrix(rnorm(m*3), ncol=3))
# Y = dist(matrix(rnorm(n*3), ncol=3))
# gwxy = gwD(X, Y)
