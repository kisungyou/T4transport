#' Barycenter of Images
#'
#' Using exact balanced optimal transport as a subroutine,
#' \code{imagebary} computes an unregularized 2-Wasserstein barycenter image
#' \eqn{X^\star} from multiple input images \eqn{X_1,\ldots,X_N}. 
#' Unlike the other image barycenter routines, this function does not use 
#' entropic regularization. Instead, it solves the barycenter problem with a 
#' robust first-order method based on mirror descent on the probability simplex.
#'
#' The algorithm treats each image as a discrete probability distribution on a
#' common \eqn{(m\times n)} grid. At each iteration, it computes exact OT dual
#' potentials \eqn{u_i} between the current barycenter iterate and each input
#' image via \code{util_dual_emd_C}. These dual potentials form a valid subgradient
#' of the barycenter objective, and a KL-mirror descent step produces a strictly
#' positive update of the barycenter weights. For numerical stability, the
#' implementation includes (i) centering of dual potentials (shift invariance),
#' (ii) gradient clipping, (iii) log-domain normalization, and (iv) optional
#' smoothing/backtracking safeguards to avoid infeasible OT calls.
#'
#' @param images a length-\eqn{N} list of same-size image matrices of size \eqn{(m\times n)}.
#' @param p an exponent for the order of the distance (default: 2). Currently, only \code{p=2}
#'   is supported (squared ground distance cost).
#' @param weights a weight of each image; if \code{NULL} (default), uniform weight is set.
#'   Otherwise, it should be a length-\eqn{N} vector of nonnegative weights.
#' @param C an optional \eqn{(mn\times mn)} ground cost matrix. If \code{NULL} (default), the
#'   squared Euclidean grid cost is used. Providing \code{C} allows using alternative ground
#'   costs (e.g., geodesic distances on a manifold discretization).
#' @param ... extra parameters including \describe{
#' \item{abstol}{stopping criterion based on \eqn{\ell_2} change of iterates (default: \code{1e-7}).}
#' \item{init.image}{an initial barycenter image (default: arithmetic mean of normalized inputs).}
#' \item{maxiter}{maximum number of mirror descent iterations (default: \code{200}).}
#' \item{step0}{initial stepsize for mirror descent (default: \code{0.5}).}
#' \item{stepschedule}{stepsize schedule; \code{"sqrt"} uses \eqn{\eta_t=\text{step0}/\sqrt{t}},
#'   and \code{"const"} uses \eqn{\eta_t=\text{step0}} (default: \code{"sqrt"}).}
#' \item{eps}{positivity floor for the barycenter and inputs; values are truncated below \code{eps}
#'   and renormalized (default: \code{1e-15}). Larger values can improve robustness.}
#' \item{smooth}{optional mixing weight toward uniform distribution after each update, used to
#'   prevent near-zero support that may cause OT infeasibility (default: \code{1e-12}). Set to
#'   \code{0} to disable.}
#' \item{clip}{\eqn{\ell_\infty} clipping threshold for the subgradient to stabilize exponentials
#'   in the KL update (default: \code{50}). Set to \code{Inf} to disable.}
#' \item{max_backtrack}{maximum number of backtracking halvings of the stepsize when an OT probe
#'   fails at the proposed update (default: \code{8}).}
#' \item{print.progress}{a logical to show iteration diagnostics (default: \code{FALSE}).}
#' }
#'
#' @return an \eqn{(m\times n)} matrix of the barycentric image.
#'
#' @examples
#' \dontrun{
#' #----------------------------------------------------------------------
#' #                       MNIST Data with Digit 3
#' #
#' # small example to compare the un- and regularized problem solutions
#' # choose only 10 images and run for 20 iterations with default penalties
#' #----------------------------------------------------------------------
#' # LOAD DATA
#' set.seed(11)
#' data(digit3)
#' dat_small = digit3[sample(1:2000, 10)]
#' 
#' # RUN
#' run_exact = imagebary(dat_small, maxiter=20)
#' run_reg14 = imagebary14C(dat_small, maxiter=20)
#' run_reg15 = imagebary15B(dat_small, maxiter=20)
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(run_exact, axes=FALSE, main="Unregularized")
#' image(run_reg14, axes=FALSE, main="Cuturi & Doucet (2014)")
#' image(run_reg15, axes=FALSE, main="Benamou et al. (2015)")
#' par(opar)
#' }
#' 
#' @concept image
#' @export
imagebary <- function(images, p = 2, weights = NULL, C = NULL, ...) {
  name.f <- "imagebary"
  params <- list(...)
  pnames <- names(params)
  
  maxiter <- as.integer(ifelse("maxiter" %in% pnames, params$maxiter, 200L))
  abstol  <- as.double(ifelse("abstol"  %in% pnames, params$abstol,  1e-7))
  step0   <- as.double(ifelse("step0"   %in% pnames, params$step0,   0.5))
  schedule <- ifelse("stepschedule" %in% pnames, params$stepschedule, "sqrt")
  eps     <- as.double(ifelse("eps" %in% pnames, params$eps, 1e-15))
  smooth  <- as.double(ifelse("smooth" %in% pnames, params$smooth, 1e-12))
  clip    <- as.double(ifelse("clip" %in% pnames, params$clip, 50))
  max_bt  <- as.integer(ifelse("max_backtrack" %in% pnames, params$max_backtrack, 8L))
  show    <- as.logical(ifelse("print.progress" %in% pnames, params$print.progress, FALSE))
  
  schedule <- match.arg(schedule, c("sqrt", "const"))
  if (p != 2) stop("* imagebary: currently only p=2 is supported.")
  if (!is.list(images) || length(images) < 1L) stop("* imagebary: 'images' must be a non-empty list.")
  if (!all(vapply(images, is.matrix, logical(1)))) stop("* imagebary: each element must be a matrix.")
  d0 <- dim(images[[1]])
  if (!all(vapply(images, function(z) all(dim(z) == d0), logical(1)))) stop("* imagebary: all images must have same size.")
  
  N <- length(images)
  w <- valid_multiple_weight(weights, N, name.f)
  w <- w / sum(w)
  
  m <- d0[1]; n <- d0[2]
  K <- m * n
  
  # cost matrix
  if (is.null(C)) {
    coordx <- seq(0, 1, length.out = n)
    coordy <- seq(1, 0, length.out = m)
    coords <- expand.grid(coordx, coordy)
    dxy <- as.matrix(stats::dist(coords))
    C <- dxy^2
  } else {
    if (!is.matrix(C) || any(dim(C) != K)) stop("* imagebary: 'C' must be (mn x mn).")
  }
  
  # normalize images onto simplex (this is essential for balanced OT)
  a_list <- lapply(images, function(X) {
    v <- as.vector(t(X))
    if (any(!is.finite(v))) stop("* imagebary: non-finite in images.")
    if (any(v < 0)) stop("* imagebary: images must be nonnegative.")
    s <- sum(v)
    if (s <= 0) stop("* imagebary: each image must have positive mass.")
    v <- v / s
    # ensure strictly positive support for stability (tiny smoothing)
    v <- pmax(v, eps)
    v / sum(v)
  })
  
  # init
  if ("init.image" %in% pnames) {
    b <- as.vector(t(params$init.image))
    if (length(b) != K || any(!is.finite(b)) || any(b < 0)) stop("* imagebary: bad init.image.")
    b <- b / sum(b)
  } else {
    b <- Reduce(`+`, a_list) / N
  }
  b <- pmax(b, eps); b <- b / sum(b)
  
  # helper: safe log-sum-exp normalization
  normalize_from_log <- function(logx) {
    mx <- max(logx)
    y <- exp(logx - mx)
    y <- y / sum(y)
    # floor + renorm
    y <- pmax(y, eps)
    y / sum(y)
  }
  
  # main loop
  for (it in seq_len(maxiter)) {
    # compute subgradient g = sum_i w_i u_i
    g <- numeric(K)
    
    ok <- TRUE
    for (i in seq_len(N)) {
      # Use tryCatch: solver can fail if b gets too sparse / tiny numerical mismatch
      du <- tryCatch(
        util_dual_emd_C(b, a_list[[i]], C, return_plan = FALSE),
        error = function(e) e
      )
      if (inherits(du, "error")) {
        ok <- FALSE
        break
      }
      ui <- as.numeric(du$u)
      
      # Center potential to remove shift ambiguity: E_b[ui]=0
      ui <- ui - sum(ui * b)
      
      g <- g + w[i] * ui
    }
    
    # choose step
    eta_base <- if (schedule == "sqrt") step0 / sqrt(it) else step0
    
    # if OT solve failed, damp and smooth b (fallback)
    if (!ok) {
      if (show) message(sprintf("* imagebary: OT failed at iter %d; smoothing/damping.", it))
      # push b slightly toward uniform and continue
      b <- (1 - 1e-2) * b + 1e-2 * (rep(1 / K, K))
      b <- pmax(b, eps); b <- b / sum(b)
      next
    }
    
    # clip gradient in linf to prevent exp overflow/underflow
    if (is.finite(clip) && clip > 0) {
      g <- pmin(pmax(g, -clip), clip)
    }
    
    # backtracking on eta if update causes solver failure next iter
    b_old <- b
    eta <- eta_base
    success <- FALSE
    
    for (bt in 0:max_bt) {
      # log-domain multiplicative update: log b_new ∝ log b - eta*g
      logb_new <- log(b_old) - eta * g
      
      b_new <- normalize_from_log(logb_new)
      
      # optional smoothing to avoid near-zero support that upsets compression
      # smooth is a tiny mixture with uniform
      if (smooth > 0) {
        b_new <- (1 - smooth) * b_new + smooth * (rep(1 / K, K))
        b_new <- pmax(b_new, eps); b_new <- b_new / sum(b_new)
      }
      
      # quick feasibility probe: try one OT call (cheap relative to full N loop)
      probe <- tryCatch(
        util_dual_emd_C(b_new, a_list[[1]], C, return_plan = FALSE),
        error = function(e) e
      )
      if (!inherits(probe, "error")) {
        b <- b_new
        success <- TRUE
        break
      }
      
      # reduce eta and try again
      eta <- eta / 2
    }
    
    if (!success) {
      # last resort: heavily smooth and continue
      if (show) message(sprintf("* imagebary: backtracking failed at iter %d; heavy smoothing.", it))
      b <- (1 - 5e-2) * b_old + 5e-2 * (rep(1 / K, K))
      b <- pmax(b, eps); b <- b / sum(b)
      next
    }
    
    # stopping
    d <- sqrt(sum((b - b_old)^2))
    if (show && (it == 1L || it %% 10L == 0L)) {
      message(sprintf("* imagebary: iter %d, eta=%.3e, ||db||_2=%.3e", it, eta, d))
    }
    if (d < abstol) break
  }
  
  return(matrix(b, m, n, byrow = TRUE))
}