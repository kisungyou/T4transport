# Barycenter of Images

Using exact balanced optimal transport as a subroutine, `imagebary`
computes an unregularized 2-Wasserstein barycenter image \\X^\star\\
from multiple input images \\X_1,\ldots,X_N\\. Unlike the other image
barycenter routines, this function does not use entropic regularization.
Instead, it solves the barycenter problem with a robust first-order
method based on mirror descent on the probability simplex.

## Usage

``` r
imagebary(images, p = 2, weights = NULL, C = NULL, ...)
```

## Arguments

- images:

  a length-\\N\\ list of same-size image matrices of size \\(m\times
  n)\\.

- p:

  an exponent for the order of the distance (default: 2). Currently,
  only `p=2` is supported (squared ground distance cost).

- weights:

  a weight of each image; if `NULL` (default), uniform weight is set.
  Otherwise, it should be a length-\\N\\ vector of nonnegative weights.

- C:

  an optional \\(mn\times mn)\\ ground cost matrix. If `NULL` (default),
  the squared Euclidean grid cost is used. Providing `C` allows using
  alternative ground costs (e.g., geodesic distances on a manifold
  discretization).

- ...:

  extra parameters including

  abstol

  :   stopping criterion based on \\\ell_2\\ change of iterates
      (default: `1e-7`).

  init.image

  :   an initial barycenter image (default: arithmetic mean of
      normalized inputs).

  maxiter

  :   maximum number of mirror descent iterations (default: `200`).

  step0

  :   initial stepsize for mirror descent (default: `0.5`).

  stepschedule

  :   stepsize schedule; `"sqrt"` uses \\\eta_t=\text{step0}/\sqrt{t}\\,
      and `"const"` uses \\\eta_t=\text{step0}\\ (default: `"sqrt"`).

  eps

  :   positivity floor for the barycenter and inputs; values are
      truncated below `eps` and renormalized (default: `1e-15`). Larger
      values can improve robustness.

  smooth

  :   optional mixing weight toward uniform distribution after each
      update, used to prevent near-zero support that may cause OT
      infeasibility (default: `1e-12`). Set to `0` to disable.

  clip

  :   \\\ell\_\infty\\ clipping threshold for the subgradient to
      stabilize exponentials in the KL update (default: `50`). Set to
      `Inf` to disable.

  max_backtrack

  :   maximum number of backtracking halvings of the stepsize when an OT
      probe fails at the proposed update (default: `8`).

  print.progress

  :   a logical to show iteration diagnostics (default: `FALSE`).

## Value

an \\(m\times n)\\ matrix of the barycentric image.

## Details

The algorithm treats each image as a discrete probability distribution
on a common \\(m\times n)\\ grid. At each iteration, it computes exact
OT dual potentials \\u_i\\ between the current barycenter iterate and
each input image via `util_dual_emd_C`. These dual potentials form a
valid subgradient of the barycenter objective, and a KL-mirror descent
step produces a strictly positive update of the barycenter weights. For
numerical stability, the implementation includes (i) centering of dual
potentials (shift invariance), (ii) gradient clipping, (iii) log-domain
normalization, and (iv) optional smoothing/backtracking safeguards to
avoid infeasible OT calls.

## Examples

``` r
if (FALSE) { # \dontrun{
#----------------------------------------------------------------------
#                       MNIST Data with Digit 3
#
# small example to compare the un- and regularized problem solutions
# choose only 10 images and run for 20 iterations with default penalties
#----------------------------------------------------------------------
# LOAD DATA
set.seed(11)
data(digit3)
dat_small = digit3[sample(1:2000, 10)]

# RUN
run_exact = imagebary(dat_small, maxiter=20)
run_reg14 = imagebary14C(dat_small, maxiter=20)
run_reg15 = imagebary15B(dat_small, maxiter=20)

# VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,3), pty="s")
image(run_exact, axes=FALSE, main="Unregularized")
image(run_reg14, axes=FALSE, main="Cuturi & Doucet (2014)")
image(run_reg15, axes=FALSE, main="Benamou et al. (2015)")
par(opar)
} # }
```
