# Gromov-Wasserstein Barycenter

Computes the Gromov–Wasserstein (GW) barycenter of a collection of
metric measure spaces. Given a list of distance matrices
\\{D^{(k)}}{k=1}^K\\ and their corresponding marginal distributions, the
function estimates a synthetic metric space whose intrinsic geometry
best represents the input collection under the GW criterion.

The GW barycenter is defined as the minimizer of a multi-measure
Gromov–Wasserstein objective, where each dataset contributes according
to a user-specified barycentric weight. Since the problem is jointly
non-convex in the barycenter metric and the coupling matrices, the
algorithm proceeds through an outer–inner iterative procedure.

## Usage

``` r
gwbary(distances, marginals = NULL, weights = NULL, num_support = 100, ...)
```

## Arguments

- distances:

  a length-\\K\\ list where each element is either an \\(N_k \times
  N_k)\\ distance matrix or an object of class `dist` representing the
  pairwise distances for each empirical measure.

- marginals:

  marginal distributions for empirical measures; if `NULL` (default),
  uniform weights are set for all measures. Otherwise, it should be a
  length-\\K\\ list where each element is a length-\\N_k\\ vector of
  nonnegative weights that sum to 1.

- weights:

  weights for each individual measure; if `NULL` (default), each measure
  is considered equally. Otherwise, it should be a length-\\K\\ vector.

- num_support:

  the number of support points \\M\\ for the barycenter (default: 100).

- ...:

  extra parameters including

  maxiter

  :   maximum number of iterations (default: 10).

  abstol

  :   stopping criterion for iterations (default: 1e-6).

  method

  :   optimization method to use; can be one of `"mm"`, `"pg"`, or
      `"fw"` (default).

## Value

A named list containing

- dist:

  an object of class `dist` representing the GW barycenter.

- weight:

  a length-\\M\\ vector of barycenter weights with all entries being
  \\1/M\\.

## Examples

``` r
if (FALSE) { # \dontrun{
#-------------------------------------------------------------------
#                           Description
#
# GW barycenter computation is quite expensive. In this example, 
# we draw a small set of empirical measures from the digit '3'
# images and compute their GW barycenter with a small number of 
# support points. The attained barycenter distance matrix is then 
# passed onto the classical MDS algorithm for visualization.
#-------------------------------------------------------------------
## GENERATE DATA
data(digits)
data_D = vector("list", length=5)
data_W = vector("list", length=5)
for (i in 1:5){
  img_now = img2measure(digits3[[i]])
  data_D[[i]] = stats::dist(img_now$support)
  data_W[[i]] = as.vector(img_now$weight)
}

## COMPUTE
bary_dist <- gwbary(data_D, marginals=data_W, num_support=100)
bary_cmd2 <- stats::cmdscale(bary_dist$dist, k=2)

## VISUALIZE
opar <- par(no.readonly=TRUE)
par(pty="s")
plot(bary_cmd2, main="GW Barycenter Embedding",
     xaxt="n", yaxt="n", pch=19, xlab="", ylab="")
par(opar)
} # }
```
