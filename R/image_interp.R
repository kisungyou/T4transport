#' Interpolation between Images 
#'
#' Given two grayscale images represented as numeric matrices of identical size,
#' compute interpolated images along a 2-Wasserstein geodesic connecting them.
#' The function interprets each image as a discrete probability distribution on
#' a common \eqn{(m\times n)} grid, computes an exact optimal transport plan, 
#' and constructs intermediate measures by pushing the
#' plan through the linear interpolation map \eqn{z=(1-t)x+t y} (displacement
#' interpolation / McCann's interpolation).
#'
#' Because the interpolated support locations generally do not coincide with
#' the original grid points, the resulting distribution is projected back onto
#' the grid by depositing transported mass to the nearest grid location. 
#' This is a simple and robust "re-binning" step, analogous in spirit to how 
#' \code{histinterp} re-bins interpolated quantile samples.
#'
#' @param image1 a grayscale image matrix of size \eqn{(m\times n)} with nonnegative entries.
#' @param image2 another grayscale image matrix of size \eqn{(m\times n)} with nonnegative entries.
#' @param t a scalar or numeric vector in \eqn{[0,1]} specifying interpolation times.
#'   \code{t=0} returns \code{image1}, \code{t=1} returns \code{image2}.
#' @param ... extra parameters including \describe{
#' \item{eps}{positivity floor applied after normalization (default: \code{1e-15}).
#'   Larger values can improve robustness.}
#' \item{abstol}{tolerance used for internal mass checks (default: \code{1e-12}).}
#' \item{print.progress}{logical; if \code{TRUE}, print basic diagnostics (default: \code{FALSE}).}
#' }
#'
#' @return
#' If \code{length(t)==1}, a single \eqn{(m\times n)} matrix representing the interpolated image.
#' If \code{length(t)>1}, a length-\code{length(t)} list of \eqn{(m\times n)} matrices.
#'
#' @examples
#' \donttest{
#' #----------------------------------------------------------------------
#' #                  Digit Interpolation between 1 and 8
#' #----------------------------------------------------------------------
#' # LOAD DATA
#' set.seed(11)
#' data(digits)
#' x1 <- digits$image[[sample(which(digits$label==1),1)]] 
#' x2 <- digits$image[[sample(which(digits$label==8),1)]]
#'
#' # COMPUTE
#' tvec <- seq(0, 1, length.out=10)
#' path <- imageinterp(x1, x2, t = tvec)
#'
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,5), pty="s")
#' for (k in 1:10){
#'   image(path[[k]], axes=FALSE, main=sprintf("t=%.2f", tvec[k]))
#' } 
#' par(opar)
#' }
#'
#' @concept image
#' @export
imageinterp <- function(image1, image2, t = 0.5, ...) {
  name.f <- "imageinterp"
  params <- list(...)
  pnames <- names(params)
  
  project <- "nearest"
  normalize <- TRUE
  C <- NULL
  
  eps <- as.double(ifelse("eps" %in% pnames, params$eps, 1e-15))
  tol <- as.double(ifelse("abstol" %in% pnames, params$abstol, 1e-12))
  project <- ifelse("project" %in% pnames, params$project, "nearest")
  show <- as.logical(ifelse("print.progress" %in% pnames, params$print.progress, FALSE))
  
  project <- match.arg(project, c("nearest"))
  
  # checks
  if (!is.matrix(image1) || !is.matrix(image2)) stop("* imageinterp: inputs must be matrices.")
  if (any(dim(image1) != dim(image2))) stop("* imageinterp: images must have the same size.")
  if (any(!is.finite(image1)) || any(!is.finite(image2))) stop("* imageinterp: non-finite values in inputs.")
  if (any(image1 < 0) || any(image2 < 0)) stop("* imageinterp: images must be nonnegative.")
  
  t <- as.numeric(t)
  if (any(t < 0 | t > 1)) stop("* imageinterp: all 't' values must lie in [0,1].")
  
  m <- nrow(image1); n <- ncol(image1)
  K <- m * n
  
  # grid coordinates (consistent with your other image routines)
  coordx <- seq(0, 1, length.out = n)
  coordy <- seq(1, 0, length.out = m)
  coords <- as.matrix(expand.grid(coordx, coordy))  # (K x 2), rows match vectorization below
  
  # vectorize in your convention
  a <- as.vector(t(image1))
  b <- as.vector(t(image2))
  
  if (normalize) {
    sa <- sum(a); sb <- sum(b)
    if (sa <= 0 || sb <= 0) stop("* imageinterp: each image must have positive mass when normalize=TRUE.")
    a <- a / sa
    b <- b / sb
    a <- pmax(a, eps); a <- a / sum(a)
    b <- pmax(b, eps); b <- b / sum(b)
  } else {
    sa <- sum(a); sb <- sum(b)
    if (abs(sa - sb) > tol) stop("* imageinterp: sums must match when normalize=FALSE (balanced OT).")
    # still guard against exact zeros if requested
    if (eps > 0) {
      a <- pmax(a, eps); a <- a / sum(a)
      b <- pmax(b, eps); b <- b / sum(b)
    }
  }
  
  # cost matrix (squared Euclidean)
  if (is.null(C)) {
    dxy <- as.matrix(stats::dist(coords))
    C <- dxy^2
  } else {
    if (!is.matrix(C) || any(dim(C) != K)) stop("* imageinterp: 'C' must be (mn x mn).")
    if (any(!is.finite(C)) || any(C < 0)) stop("* imageinterp: 'C' must be finite and nonnegative.")
  }
  
  # OT plan between endpoints (exact)
  if (show) message("* imageinterp: solving exact OT plan...")
  G <- util_plan_emd_C(a, b, C)
  
  # indices of nonzero flows (sparse traversal)
  nz <- which(G > 0, arr.ind = TRUE)
  if (nrow(nz) == 0L) {
    # Degenerate: return something sensible
    if (length(t) == 1L) return(matrix(a, m, n, byrow = TRUE))
    out <- vector("list", length(t))
    for (k in seq_along(t)) out[[k]] <- matrix((1 - t[k]) * a + t[k] * b, m, n, byrow = TRUE)
    return(out)
  }
  
  ii <- nz[, 1]
  jj <- nz[, 2]
  flow <- G[nz]
  
  xi <- coords[ii, , drop = FALSE]  # (nnz x 2)
  yj <- coords[jj, , drop = FALSE]  # (nnz x 2)
  
  # helper: project mass at (x,y) to nearest grid index
  # Since the grid is regular, nearest index can be computed by rounding.
  # coords come from expand.grid(coordx, coordy):
  #   x-index corresponds to coordx (length n)
  #   y-index corresponds to coordy (length m)
  nearest_grid_index <- function(zx, zy) {
    # map zx to column index in 1..n
    jx <- round((zx - coordx[1]) / (coordx[n] - coordx[1]) * (n - 1)) + 1
    jx <- pmin(pmax(jx, 1), n)
    
    # map zy to row index in 1..m
    iy <- round((coordy[1] - zy) / (coordy[1] - coordy[m]) * (m - 1)) + 1
    iy <- pmin(pmax(iy, 1), m)
    
    # convert (row, col) to vector index consistent with as.vector(t(.))
    # as.vector(t(M)) stacks rows of M (row-major), so index = (row-1)*n + col
    (iy - 1) * n + jx
  }
  
  make_image_at_t <- function(tau) {
    if (tau <= 0) return(matrix(a, m, n, byrow = TRUE))
    if (tau >= 1) return(matrix(b, m, n, byrow = TRUE))
    
    # interpolated positions for each transported mass atom
    z <- (1 - tau) * xi + tau * yj
    idx <- nearest_grid_index(z[, 1], z[, 2])
    
    out <- numeric(K)
    
    # robust 1-way accumulation: out[k] = sum_{ell: idx[ell]=k} flow[ell]
    acc <- rowsum(flow, idx, reorder = FALSE)  # matrix
    out[as.integer(rownames(acc))] <- acc[, 1]
    
    # numerical cleanup / normalization back to simplex
    out <- pmax(out, 0)
    s <- sum(out)
    if (!is.finite(s) || s <= 0) {
      out <- (1 - tau) * a + tau * b
      out <- pmax(out, eps); out <- out / sum(out)
    } else {
      out <- out / s
      out <- pmax(out, eps); out <- out / sum(out)
    }
    
    matrix(out, m, n, byrow = TRUE)
  }
  
  
  if (length(t) == 1L) {
    return(make_image_at_t(t))
  } else {
    out_list <- vector("list", length(t))
    for (k in seq_along(t)) out_list[[k]] <- make_image_at_t(t[k])
    return(out_list)
  }
}
