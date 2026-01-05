#' Wasserstein Distance between Two Images
#'
#' @description
#' Given two grayscale images represented as numeric matrices, compute their
#' Wasserstein distance using an exact balanced optimal transport solver. 
#' Each image is interpreted as a discrete probability distribution on a common \eqn{(m\times n)} grid. 
#' The ground cost is defined using the Euclidean distance between grid locations.
#'
#' @param x a grayscale image matrix of size \eqn{(m\times n)} with nonnegative entries.
#' @param y a grayscale image matrix of size \eqn{(m\times n)} with nonnegative entries.
#' @param p an exponent for the order of the distance (default: 2). 
#'
#' @return a list containing \describe{
#' \item{distance}{the Wasserstein distance \eqn{W_p(x,y)}.}
#' \item{plan}{the optimal transport plan matrix of size \eqn{(mn\times mn)}.}
#' }
#'
#' @examples
#' \donttest{
#' #----------------------------------------------------------------------
#' #                       Small MNIST-like Example
#' #----------------------------------------------------------------------
#' # DATA
#' data(digit3)
#' x <- digit3[[1]]
#' y <- digit3[[2]]
#'
#' # COMPUTE
#' W1 <- imagedist(x, y, p=1)
#' W2 <- imagedist(x, y, p=2)
#' 
#' # SHOW RESULTS
#' print(paste0("Wasserstein-1 distance: ", round(W1$distance,4)))
#' print(paste0("Wasserstein-2 distance: ", round(W2$distance,4)))
#' }
#'
#' @concept image
#' @export
imagedist <- function(x, y, p = 2) {
  name.f <- "imagedist"
  normalize = TRUE
  return_plan = TRUE
  eps <- 0
  tol <- 1e-12
  C <- NULL
  
  par_p = as.double(p)
  if (par_p < 1) stop("* imagedist: 'p' must be at least 1.")
  if (!is.matrix(x) || !is.matrix(y)) stop("* imagedist: 'x' and 'y' must be matrices.")
  if (any(dim(x) != dim(y))) stop("* imagedist: 'x' and 'y' must have the same size.")
  if (any(!is.finite(x)) || any(!is.finite(y))) stop("* imagedist: non-finite values in inputs.")
  if (any(x < 0) || any(y < 0)) stop("* imagedist: images must be nonnegative.")
  
  m <- nrow(x); n <- ncol(x)
  K <- m * n
  
  # vectorize in the same convention as your barycenter code
  a <- as.vector(t(x))
  b <- as.vector(t(y))
  
  # normalize or require equal mass
  if (normalize) {
    sa <- sum(a); sb <- sum(b)
    if (sa <= 0 || sb <= 0) stop("* imagedist: each image must have positive total mass when normalize=TRUE.")
    a <- a / sa
    b <- b / sb
    if (eps > 0) {
      a <- pmax(a, eps); a <- a / sum(a)
      b <- pmax(b, eps); b <- b / sum(b)
    }
  } else {
    sa <- sum(a); sb <- sum(b)
    if (abs(sa - sb) > tol) stop("* imagedist: sum(x) and sum(y) must match for balanced OT when normalize=FALSE.")
  }
  
  # cost matrix
  if (is.null(C)) {
    coordx <- seq(0, 1, length.out = n)
    coordy <- seq(1, 0, length.out = m)
    coords <- expand.grid(coordx, coordy)
    dxy <- as.matrix(stats::dist(coords))
    C <- dxy^par_p
  } else {
    if (!is.matrix(C) || any(dim(C) != K)) stop("* imagedist: 'C' must be (mn x mn).")
    if (any(!is.finite(C)) || any(C < 0)) stop("* imagedist: 'C' must be finite and nonnegative.")
  }
  
  if (!return_plan) {
    # Use dual call to avoid materializing plan if you later want speed.
    # But you already have plan solver; distance computed from plan is exact too.
    G <- util_plan_emd_C(a, b, C)
    cost <- sum(G * C)
    # W2 = sqrt(cost) because C is squared distances
    return(max(0, cost)^(1/p))
    return(sqrt(max(0, cost)))
  } else {
    G <- util_plan_emd_C(a, b, C)
    cost <- sum(G * C)
    return(list(distance = max(0, cost)^(1/p), plan = G))
  }
}
