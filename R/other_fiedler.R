#' Compute the fiedler vector of a point cloud
#' 
#' Given a point cloud \eqn{X \in \mathbf{R}^{N \times P}}, this function
#' constructs a fully connected weighted graph using an RBF (Gaussian) kernel
#' with bandwidth chosen by the median heuristic, forms the unnormalized graph
#' Laplacian, and returns the corresponding Fiedler vector, which is the eigenvector
#' associated to the second smallest eigenvalue of the Laplacian.
#'
#' @param X An \eqn{(N\times P)} matrix of row observations.
#' @param normalize Logical; if \code{TRUE} (default), the Fiedler vector is
#'   rescaled to lie in \eqn{[0,1]} by subtracting its minimum and dividing by
#'   its range, mimicking the normalization convention in the corresponding
#'   Python implementation. If \code{FALSE}, the raw eigenvector is returned.
#'
#' @return A numeric vector of length \eqn{N} containing the Fiedler values
#'   associated with each point in the input point cloud. If \code{normalize = TRUE},
#'   the entries are in the interval \eqn{[0,1]}.
#'   
#' @examples
#' #-------------------------------------------------------------------
#' #                           Description
#' #
#' # Use 'iris' dataset to compute fiedler vector. 
#' # The dataset is visualized in R^2 using PCA
#' #-------------------------------------------------------------------
#' # load dataset
#' X = as.matrix(iris[,1:4])
#' 
#' # PCA preprocessing
#' X2d = X%*%eigen(cov(X))$vectors[,1:2]
#' 
#' # compute fiedler vector
#' fied_vec = fiedler(X2d, normalize=TRUE)
#' 
#' # plot 
#' opar <- par(no.readonly=TRUE)
#' plot(X2d, col=rainbow(150)[as.numeric(cut(fied_vec, breaks=150))], 
#'      pch=19, xlab="PC 1", ylab="PC 2",
#'      main="Fiedler vector on Iris dataset (PCA-reduced)")
#' par(opar)
#' 
#' @concept other
#' @export
fiedler <- function(X, normalize=TRUE){
  # check the input
  if (is.vector(X)){X = matrix(X, ncol=1)}
  if (!is.matrix(X)){    stop("* fiedler : input 'X' should be a matrix.")  }
  
  # computation
  pair_D = util_pairwise_dist(X) # pairwise distance
  sigma = stats::median(as.vector(pair_D[upper.tri(pair_D)]))
  if (sigma <= .Machine$double.eps) {
    stop("* fiedler : median distance is zero; cannot define RBF kernel.")
  }
  W = exp(- (pair_D^2) / (2*sigma^2)) # affinity matrix
  
  # Degree vector and unnormalized Laplacian
  deg <- rowSums(W)
  L   <- diag(deg) - W
  
  # Eigen-decomposition (L is symmetric)
  eig <- eigen(L, symmetric = TRUE)
  
  # Sort eigenvalues ascending
  ord    <- order(eig$values)
  vals   <- eig$values[ord]
  vecs   <- eig$vectors[, ord, drop = FALSE]
  
  # Fiedler vector = eigenvector for 2nd smallest eigenvalue
  fiedler <- Re(vecs[, 2])
  
  # python convention
  if (normalize){
    # Rescale to [0, 1] as in the Python code
    f_shift <- fiedler - min(fiedler)
    v       <- f_shift / max(f_shift)
    
    return(v)
  } else {
    return(as.vector(fiedler))
  }
}
