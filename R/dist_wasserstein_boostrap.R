#' Wasserstein Distance Estimation with Boostrapping
#' 
#' This function computes the \eqn{\mathcal{W}_p} distance between two empirical measures 
#' using bootstrap in order to quantify the uncertainty of the estimation.
#' 
#' @param X an \eqn{(M\times P)} matrix of row observations.
#' @param Y an \eqn{(N\times P)} matrix of row observations.
#' @param p an exponent for the order of the distance (default: 2).
#' @param B number of bootstrap samples (default: 500).
#' @param wx a length-\eqn{M} marginal density that sums to \eqn{1}. If \code{NULL} (default), uniform weight is set.
#' @param wy a length-\eqn{N} marginal density that sums to \eqn{1}. If \code{NULL} (default), uniform weight is set.
#' 
#' @return a named list containing\describe{
#' \item{distance}{\eqn{\mathcal{W}_p} distance value.}
#' \item{boot_samples}{a length-\eqn{B} vector of bootstrap samples.}
#' }
#' 
#' @examples
#' \donttest{
#' #-------------------------------------------------------------------
#' #  Boostrapping Wasserstein Distance between Two Bivariate Normals
#' #
#' # * class 1 : samples from Gaussian with mean=(-5, 0)
#' # * class 2 : samples from Gaussian with mean=(+5, 0)
#' #-------------------------------------------------------------------
#' ## SMALL EXAMPLE
#' m = round(runif(1, min=50, max=100))
#' n = round(runif(1, min=50, max=100))
#' X = matrix(rnorm(m*2), ncol=2) # m obs. for X
#' Y = matrix(rnorm(n*2), ncol=2) # n obs. for Y
#' 
#' X[,1] = X[,1] - 5
#' Y[,1] = Y[,1] + 5
#' 
#' ## COMPUTE THE BOOTSTRAP SAMPLES
#' boots = wassboot(X, Y, B=1000)
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' hist(boots$boot_samples, xlab="Estimates", main="Bootstrap Samples")
#' abline(v=boots$distance, lwd=2, col="blue")
#' abline(v=mean(boots$boot_samples), lwd=2, col="red")
#' abline(v=10, col="cyan", lwd=2)
#' legend("topright", c("ground truth","estimate","bootstrap mean"), 
#'         col=c("cyan","blue","red"), lwd=2)
#' par(opar)
#' }
#' 
#' @concept dist
#' @export
wassboot <- function(X, Y, p=2, B=500, wx=NULL, wy=NULL){
  ## CHECK INPUTS
  if (is.vector(X)){
    X = matrix(X, ncol=1)
  }
  if (is.vector(Y)){
    Y = matrix(Y, ncol=1)
  }
  if (!is.matrix(X)){    stop("* wassboot : input 'X' should be a matrix.")  }
  if (!is.matrix(Y)){    stop("* wassboot : input 'Y' should be a matrix.")  }
  if (base::ncol(X)!=base::ncol(Y)){
    stop("* wassboot : input 'X' and 'Y' should be of same dimension.")
  }
  m = base::nrow(X)
  n = base::nrow(Y)
  
  wxname = paste0("'",deparse(substitute(wx)),"'")
  wyname = paste0("'",deparse(substitute(wy)),"'")
  fname  = "wassboot"
  
  par_wx = valid_single_marginal(wx, m, fname)
  par_wy = valid_single_marginal(wy, n, fname) #valid_weight(wy, n, wyname, fname)
  par_p  = max(1, as.double(p))
  par_D  = as.matrix(compute_pdist2(X, Y))
  par_B = max(5, round(B))
  
  # COMPUTE - ESTIMATE
  out_estimate <- wass_lp(par_D, par_p, par_wx, par_wy)$distance
  
  # COMPUTE - BOOTSTRAP SAMPLES
  out_boots = rep(0, par_B)
  for (it in seq_len(par_B)){
    out_boots[it] = aux_wasserstein_boot(par_D, par_p, par_wx, par_wy)
  }
  
  # RETURN
  return(list(distance=out_estimate, 
              boot_samples=out_boots))
}



# auxiliary function ------------------------------------------------------
#' @keywords internal
#' @noRd
aux_wasserstein_boot <- function(D, p, wx, wy){
  # normalize the weights just in case
  norm_wx = wx/base::sum(wx)
  norm_wy = wy/base::sum(wy)
  
  # empirical measure size
  nx = length(norm_wx)
  ny = length(norm_wy)
  
  # sampling
  sample_x = sample(seq_len(nx), size=nx, replace=TRUE, prob=norm_wx)
  sample_y = sample(seq_len(ny), size=ny, replace=TRUE, prob=norm_wy)
  
  # tabularize
  table_x = table(sample_x)
  table_y = table(sample_y)
  
  # from table: indices
  select_x = as.numeric(names(table_x))
  select_y = as.numeric(names(table_y))
  
  # from table: frequency
  freq_x = as.vector(table_x); freq_x = freq_x/base::sum(freq_x)
  freq_y = as.vector(table_y); freq_y = freq_y/base::sum(freq_y)
  
  # compute
  sub_D = D[select_x, select_y]
  out_estimate <- wass_lp(sub_D, p, freq_x, freq_y)$distance
  return(out_estimate)
}