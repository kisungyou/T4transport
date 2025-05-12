#' Sliced Wasserstein Distance
#' 
#' @description
#' Sliced Wasserstein (SW) Distance is a popular alternative to the standard Wasserstein distance due to its computational 
#' efficiency on top of nice theoretical properties. For the \eqn{d}-dimensional probability 
#' measures \eqn{\mu} and \eqn{\nu}, the SW distance is defined as 
#' \deqn{\mathcal{SW}_p (\mu, \nu) = 
#' \left( \int_{\mathbb{S}^{d-1}} \mathcal{W}_p^p (
#' \langle \theta, \mu\rangle, \langle \theta, \nu \rangle) d\lambda (\theta) \right)^{1/p},}
#' where \eqn{\mathbb{S}^{d-1}} is the \eqn{(d-1)}-dimensional unit hypersphere and 
#' \eqn{\lambda} is the uniform distribution on \eqn{\mathbb{S}^{d-1}}. Practically, 
#' it is computed via Monte Carlo integration.
#' 
#' @param X an \eqn{(M\times P)} matrix of row observations.
#' @param Y an \eqn{(N\times P)} matrix of row observations.
#' @param p an exponent for the order of the distance (default: 2).
#' @param ... extra parameters including \describe{
#' \item{num_proj}{the number of Monte Carlo samples for SW computation (default: 496).}
#' }
#' 
#' @return a named list containing\describe{
#' \item{distance}{\eqn{\mathcal{SW}_p} distance value.}
#' \item{projdist}{a length-\code{num_proj} vector of projected univariate distances.}
#' }
#' 
#' @examples
#' \donttest{
#' #-------------------------------------------------------------------
#' #  Sliced-Wasserstein Distance between Two Bivariate Normal
#' #
#' # * class 1 : samples from Gaussian with mean=(-1, -1)
#' # * class 2 : samples from Gaussian with mean=(+1, +1)
#' #-------------------------------------------------------------------
#' # SMALL EXAMPLE
#' set.seed(100)
#' m = 20
#' n = 30
#' X = matrix(rnorm(m*2, mean=-1),ncol=2) # m obs. for X
#' Y = matrix(rnorm(n*2, mean=+1),ncol=2) # n obs. for Y
#' 
#' # COMPUTE THE SLICED-WASSERSTEIN DISTANCE
#' outsw <- swdist(X, Y, num_proj=100)
#' 
#' # VISUALIZE
#' # prepare ingredients for plotting
#' plot_x = 1:1000
#' plot_y = base::cumsum(outsw$projdist)/plot_x
#' 
#' # draw
#' opar <- par(no.readonly=TRUE)
#' plot(plot_x, plot_y, type="b", cex=0.1, lwd=2,
#'      xlab="number of MC samples", ylab="distance",
#'      main="Effect of MC Sample Size")
#' abline(h=outsw$distance, col="red", lwd=2)
#' legend("bottomright", legend="SW Distance", 
#'        col="red", lwd=2)
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{rabin_2012_WassersteinBarycenterIts}{T4transport}
#' 
#' @concept dist
#' @name swdist
#' @rdname swdist
#' @export
swdist <- function(X, Y, p=2, ...){
  ## INPUTS : EXPLICIT
  if (is.vector(X)){
    X = matrix(X, ncol=1)
  }
  if (is.vector(Y)){
    Y = matrix(Y, ncol=1)
  }
  if (!is.matrix(X)){    stop("* swdist : input 'X' should be a matrix.")  }
  if (!is.matrix(Y)){    stop("* swdist : input 'Y' should be a matrix.")  }
  if (base::ncol(X)!=base::ncol(Y)){
    stop("* swdist : input 'X' and 'Y' should be of same dimension.")
  }
  par_p     = max(1, as.double(p))
  
  ## INPUTS : IMPLICIT
  params = list(...)
  pnames = names(params)
  
  if ("num_proj"%in%pnames){
    par_niter = max(1, round(params$num_proj))
  } else {
    par_niter = 496
  }
  
  ## COMPUTATION
  # base parameter
  m = base::nrow(X)
  n = base::nrow(Y)
  
  # case branching : univariate vs. multivariate
  if (ncol(X)==1){
    distval = swdist_univariate(as.vector(X), as.vector(Y), par_p)
    output  = list(
      distance=distval,
      projdist=NA
    )
    return(output)
  } else {
    distvec = rep(0, par_niter)
    for (it in 1:par_niter){
      # random projection
      randproj <- swdist_projection(ncol(X))
      projX <- as.vector(X%*%randproj)
      projY <- as.vector(Y%*%randproj)
      
      # computation
      distvec[it] = swdist_univariate(projX, projY, par_p)
    }
    output = list(
      distance=base::mean(distvec),
      projdist=as.vector(distvec)
    )
  }
}


# auxiliary ---------------------------------------------------------------
#' @keywords internal
#' @noRd
swdist_univariate <- function(vec1, vec2, p){
  # grid
  npts  = 1000
  eeps  = sqrt(.Machine$double.eps)
  seq_x = seq(from=eeps, to=(1-eeps), length.out=npts)
  
  # get ecdfs
  ecdf_1 = stats::ecdf(vec1)
  ecdf_2 = stats::ecdf(vec2)
  
  # quantile values
  quantile_1 = as.vector(stats::quantile(ecdf_1, seq_x))
  quantile_2 = as.vector(stats::quantile(ecdf_2, seq_x))
  
  # compute
  return(as.double(ecdf_pdist(seq_x, quantile_1, quantile_2, p)))
}

#' @keywords internal
#' @noRd
swdist_projection <- function(dim){
  randproj <- stats::rnorm(dim)
  randproj <- randproj/sqrt(sum(randproj^2))
  return(matrix(randproj, ncol=1))
}
