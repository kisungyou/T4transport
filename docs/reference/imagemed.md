# Wasserstein Median of Images

Using exact balanced optimal transport as a subroutine, `imagemed`
computes an unregularized 2-Wasserstein geometric median image
\\X^\dagger\\ from multiple input images \\X_1,\ldots,X_N\\. The
Wasserstein median is defined as a minimizer of the (weighted) sum of
Wasserstein distances, \$\$ \arg\min\_{X} \sum\_{i=1}^N w_i\\ W_2(X,
X_i). \$\$

## Usage

``` r
imagemed(images, weights = NULL, C = NULL, ...)
```

## Arguments

- images:

  a length-\\N\\ list of same-size grayscale image matrices of size
  \\(m\times n)\\.

- weights:

  a weight of each image; if `NULL` (default), uniform weight is set.
  Otherwise, it should be a length-\\N\\ vector of nonnegative weights.

- C:

  an optional \\(mn\times mn)\\ ground cost matrix (squared distances).
  If `NULL` (default), the squared Euclidean grid cost is used.

- ...:

  extra parameters including

  maxiter

  :   maximum number of IRLS outer iterations (default: `30`).

  abstol

  :   stopping tolerance based on \\\ell_2\\ change of iterates
      (default: `1e-6`).

  delta

  :   small positive number to avoid division by zero in IRLS weights
      (default: `1e-8`).

  init.image

  :   initial median iterate (default: unweighted barycenter via
      `imagebary` with a small number of iterations).

  init.bary.iter

  :   iterations for the default initialization barycenter (default:
      `10`).

  bary.maxiter

  :   maximum iterations for each barycenter subproblem (default:
      `200`).

  bary.abstol

  :   tolerance for each barycenter subproblem (default: `1e-7`).

  bary.step0

  :   initial step size for barycenter subproblem (default: `0.5`).

  bary.stepschedule

  :   `"sqrt"` or `"const"` for barycenter subproblem (default:
      `"sqrt"`).

  bary.eps

  :   positivity floor used inside barycenter (default: `1e-15`).

  bary.smooth

  :   smoothing used inside barycenter (default: `1e-12`).

  bary.clip

  :   gradient clipping used inside barycenter (default: `50`).

  bary.max_backtrack

  :   backtracking cap used inside barycenter (default: `8`).

  print.progress

  :   logical; if `TRUE`, print iteration diagnostics (default:
      `FALSE`).

## Value

an \\(m\times n)\\ matrix of the median.

## Details

Unlike Wasserstein barycenters (which minimize squared distances), the
median is a robust notion of centrality. This function solves the
problem with an iterative reweighted least squares (IRLS) scheme (a
Wasserstein analogue of Weiszfeld's algorithm). Each outer iteration
updates weights based on current distances and then solves a weighted
Wasserstein barycenter problem: \$\$ \alpha_i^{(k)} \propto
\frac{w_i}{\max(W_2(X^{(k)},X_i),\delta)}, \qquad X^{(k+1)} = \arg\min_X
\sum\_{i=1}^N \alpha_i^{(k)}\\ W_2^2(X, X_i). \$\$

The barycenter subproblem is solved by
[`imagebary`](https://www.kisungyou.com/T4transport/reference/imagebary.md)
(mirror descent with exact OT dual subgradients). Distances \\W_2\\ are
computed by exact EMD plans under the same squared ground cost.

## Examples

``` r
if (FALSE) { # \dontrun{
#----------------------------------------------------------------------
#                             MNIST Example
#
# Use 6 images from digit '8' and 4 images from digit '1'.
# The median should look closer to the shape of '8'.
#----------------------------------------------------------------------
# DATA PREP
set.seed(11)
data(digits)
dat_8 = digits$image[sample(which(digits$label==8), 6)]
dat_1 = digits$image[sample(which(digits$label==1), 4)]
dat_all = c(dat_8, dat_1)

# COMPUTE BARYCENTER AND MEDIAN
img_bary = imagebary(dat_all, maxiter=50)
img_med  = imagemed(dat_all, maxiter=50)

# VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2), pty="s")
image(img_bary, axes=FALSE, main="Barycenter")
image(img_med,  axes=FALSE, main="Median")
par(opar)
} # }
```
