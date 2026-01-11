#' Wasserstein Median of Images
#'
#' Using exact balanced optimal transport as a subroutine, \code{imagemed}
#' computes an unregularized 2-Wasserstein geometric median image \eqn{X^\dagger}
#' from multiple input images \eqn{X_1,\ldots,X_N}. The Wasserstein median is
#' defined as a minimizer of the (weighted) sum of Wasserstein distances,
#' \deqn{ \arg\min_{X} \sum_{i=1}^N w_i\, W_2(X, X_i). }
#'
#' Unlike Wasserstein barycenters (which minimize squared distances), the median
#' is a robust notion of centrality. This function solves the problem with an
#' iterative reweighted least squares (IRLS) scheme (a Wasserstein analogue of
#' Weiszfeld's algorithm). Each outer iteration updates weights based on current
#' distances and then solves a weighted Wasserstein barycenter problem:
#' \deqn{ \alpha_i^{(k)} \propto \frac{w_i}{\max(W_2(X^{(k)},X_i),\delta)}, \qquad
#'        X^{(k+1)} = \arg\min_X \sum_{i=1}^N \alpha_i^{(k)}\, W_2^2(X, X_i). }
#'
#' The barycenter subproblem is solved by \code{\link{imagebary}} (mirror descent
#' with exact OT dual subgradients). Distances \eqn{W_2} are computed by exact
#' EMD plans under the same squared ground cost.
#'
#' @param images a length-\eqn{N} list of same-size grayscale image matrices of size \eqn{(m\times n)}.
#' @param weights a weight of each image; if \code{NULL} (default), uniform weight is set.
#'   Otherwise, it should be a length-\eqn{N} vector of nonnegative weights.
#' @param C an optional \eqn{(mn\times mn)} ground cost matrix (squared distances). If \code{NULL}
#'   (default), the squared Euclidean grid cost is used.
#' @param ... extra parameters including \describe{
#' \item{maxiter}{maximum number of IRLS outer iterations (default: \code{30}).}
#' \item{abstol}{stopping tolerance based on \eqn{\ell_2} change of iterates (default: \code{1e-6}).}
#' \item{delta}{small positive number to avoid division by zero in IRLS weights
#'   (default: \code{1e-8}).}
#' \item{init.image}{initial median iterate (default: unweighted barycenter via \code{imagebary}
#'   with a small number of iterations).}
#' \item{init.bary.iter}{iterations for the default initialization barycenter (default: \code{10}).}
#' \item{bary.maxiter}{maximum iterations for each barycenter subproblem (default: \code{200}).}
#' \item{bary.abstol}{tolerance for each barycenter subproblem (default: \code{1e-7}).}
#' \item{bary.step0}{initial step size for barycenter subproblem (default: \code{0.5}).}
#' \item{bary.stepschedule}{\code{"sqrt"} or \code{"const"} for barycenter subproblem (default: \code{"sqrt"}).}
#' \item{bary.eps}{positivity floor used inside barycenter (default: \code{1e-15}).}
#' \item{bary.smooth}{smoothing used inside barycenter (default: \code{1e-12}).}
#' \item{bary.clip}{gradient clipping used inside barycenter (default: \code{50}).}
#' \item{bary.max_backtrack}{backtracking cap used inside barycenter (default: \code{8}).}
#' \item{print.progress}{logical; if \code{TRUE}, print iteration diagnostics (default: \code{FALSE}).}
#' }
#'
#' @return an \eqn{(m\times n)} matrix of the median.
#'
#' @examples
#' \dontrun{
#' #----------------------------------------------------------------------
#' #                             MNIST Example
#' #
#' # Use 6 images from digit '8' and 4 images from digit '1'.
#' # The median should look closer to the shape of '8'.
#' #----------------------------------------------------------------------
#' # DATA PREP
#' set.seed(11)
#' data(digits)
#' dat_8 = digits$image[sample(which(digits$label==8), 6)]
#' dat_1 = digits$image[sample(which(digits$label==1), 4)]
#' dat_all = c(dat_8, dat_1)
#' 
#' # COMPUTE BARYCENTER AND MEDIAN
#' img_bary = imagebary(dat_all, maxiter=50)
#' img_med  = imagemed(dat_all, maxiter=50)
#'
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(img_bary, axes=FALSE, main="Barycenter")
#' image(img_med,  axes=FALSE, main="Median")
#' par(opar)
#' }
#'
#' @concept image
#' @export
imagemed <- function(images, weights = NULL, C = NULL, ...) {
  name.f <- "imagemed"
  params <- list(...)
  pnames <- names(params)
  
  # --- IRLS controls
  maxiter <- as.integer(ifelse("maxiter" %in% pnames, params$maxiter, 30L))
  abstol  <- as.double(ifelse("abstol"  %in% pnames, params$abstol,  1e-6))
  delta   <- as.double(ifelse("delta"   %in% pnames, params$delta,   1e-8))
  show    <- as.logical(ifelse("print.progress" %in% pnames, params$print.progress, FALSE))
  retw    <- FALSE
  
  init_bary_iter <- as.integer(ifelse("init.bary.iter" %in% pnames, params$init.bary.iter, 10L))
  
  # --- barycenter subproblem controls (passed explicitly)
  bary_maxiter <- as.integer(ifelse("bary.maxiter" %in% pnames, params$bary.maxiter, 200L))
  bary_abstol  <- as.double(ifelse("bary.abstol"  %in% pnames, params$bary.abstol,  1e-7))
  bary_step0   <- as.double(ifelse("bary.step0"   %in% pnames, params$bary.step0,   0.5))
  bary_sched   <- ifelse("bary.stepschedule" %in% pnames, params$bary.stepschedule, "sqrt")
  bary_eps     <- as.double(ifelse("bary.eps" %in% pnames, params$bary.eps, 1e-15))
  bary_smooth  <- as.double(ifelse("bary.smooth" %in% pnames, params$bary.smooth, 1e-12))
  bary_clip    <- as.double(ifelse("bary.clip" %in% pnames, params$bary.clip, 50))
  bary_bt      <- as.integer(ifelse("bary.max_backtrack" %in% pnames, params$bary.max_backtrack, 8L))
  
  bary_sched <- match.arg(bary_sched, c("sqrt", "const"))
  
  # --- basic checks
  if (!is.list(images) || length(images) < 1L) stop("* imagemed: 'images' must be a non-empty list.")
  if (!all(vapply(images, is.matrix, logical(1)))) stop("* imagemed: each element must be a matrix.")
  d0 <- dim(images[[1]])
  if (!all(vapply(images, function(z) all(dim(z) == d0), logical(1)))) stop("* imagemed: all images must have identical dimensions.")
  
  N <- length(images)
  w0 <- valid_multiple_weight(weights, N, name.f)
  w0 <- w0 / sum(w0)
  
  m <- d0[1]; n <- d0[2]
  K <- m * n
  
  # --- cost matrix (squared Euclidean on grid) if not provided
  if (is.null(C)) {
    coordx <- seq(0, 1, length.out = n)
    coordy <- seq(1, 0, length.out = m)
    coords <- expand.grid(coordx, coordy)
    dxy <- as.matrix(stats::dist(coords))
    C <- dxy^2
  } else {
    if (!is.matrix(C) || any(dim(C) != K)) stop("* imagemed: 'C' must be (mn x mn).")
  }
  
  # --- normalize all images to simplex vectors once (for distance evaluation)
  a_list <- lapply(images, function(X) {
    v <- as.vector(t(X))
    if (any(!is.finite(v))) stop("* imagemed: non-finite in images.")
    if (any(v < 0)) stop("* imagemed: images must be nonnegative.")
    s <- sum(v)
    if (s <= 0) stop("* imagemed: each image must have positive total mass.")
    v <- v / s
    v <- pmax(v, bary_eps)
    v / sum(v)
  })
  
  # --- initialize
  if ("init.image" %in% pnames) {
    b <- as.vector(t(params$init.image))
    if (length(b) != K || any(!is.finite(b)) || any(b < 0)) stop("* imagemed: bad init.image.")
    b <- b / sum(b)
    b <- pmax(b, bary_eps); b <- b / sum(b)
  } else {
    # a cheap, stable initialization: a few iterations of unweighted barycenter
    b_init_mat <- imagebary(images,
                            p = 2, weights = w0, C = C,
                            maxiter = init_bary_iter,
                            abstol  = bary_abstol,
                            step0   = bary_step0,
                            stepschedule = bary_sched,
                            eps = bary_eps, smooth = bary_smooth,
                            clip = bary_clip, max_backtrack = bary_bt,
                            print.progress = FALSE)
    b <- as.vector(t(b_init_mat))
    b <- b / sum(b)
    b <- pmax(b, bary_eps); b <- b / sum(b)
  }
  
  # helper: compute W2 distance between current b and a_i using exact plan
  # W2^2 = <G, C>, W2 = sqrt(W2^2)
  w2dist_to <- function(a) {
    G <- util_plan_emd_C(b, a, C)
    cost <- sum(G * C)
    sqrt(max(0, cost))
  }
  
  # IRLS outer loop
  alpha <- w0
  for (it in seq_len(maxiter)) {
    b_old <- b
    
    # 1) distances r_i = W2(b, a_i)
    r <- numeric(N)
    for (i in seq_len(N)) {
      # tryCatch to avoid hard failure; if fail, inflate r_i to reduce its weight
      ri <- tryCatch(w2dist_to(a_list[[i]]), error = function(e) NA_real_)
      if (!is.finite(ri)) ri <- Inf
      r[i] <- ri
    }
    
    # 2) IRLS weights alpha_i ∝ w_i / max(r_i, delta)
    denom <- pmax(r, delta)
    alpha <- w0 / denom
    alpha <- pmax(alpha, 0)
    if (sum(alpha) <= 0) {
      alpha <- w0
    } else {
      alpha <- alpha / sum(alpha)
    }
    
    # 3) barycenter subproblem: minimize sum_i alpha_i W2^2(., a_i)
    b_mat <- imagebary(images,
                       p = 2,
                       weights = alpha,
                       C = C,
                       maxiter = bary_maxiter,
                       abstol  = bary_abstol,
                       step0   = bary_step0,
                       stepschedule = bary_sched,
                       eps = bary_eps,
                       smooth = bary_smooth,
                       clip = bary_clip,
                       max_backtrack = bary_bt,
                       print.progress = FALSE)
    
    b <- as.vector(t(b_mat))
    b <- b / sum(b)
    b <- pmax(b, bary_eps); b <- b / sum(b)
    
    # stopping
    d <- sqrt(sum((b - b_old)^2))
    
    if (show) {
      # objective proxy: sum w_i * W2(b, a_i)
      obj <- sum(w0 * pmax(r, 0))
      message(sprintf("* imagemed: iter %d, ||db||_2=%.3e, proxy obj=%.6f", it, d, obj))
    }
    
    if (d < abstol) break
  }
  
  out_mat <- matrix(b, m, n, byrow = TRUE)
  if (!retw) return(out_mat)
  return(list(median = out_mat, weights = alpha))
}