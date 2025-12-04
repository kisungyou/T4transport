# Free-Support Median by Weiszfeld Update with Barycentric Projection

For a collection of empirical measures \\\lbrace
\mu_k\rbrace\_{k=1}^K\\, the free-support Wasserstein median, a
minimizer to the following functional \$\$ \mathcal{F}(\nu) =
\sum\_{k=1}^K w_k \mathcal{W}\_2 (\nu, \mu_k ), \$\$ is computed using
the OT-adapted version of the Weiszfeld algorithm using the barycentric
projection as a means to recover an optimal displacement map.

## Usage

``` r
rmedWB(atoms, marginals = NULL, weights = NULL, num_support = 100, ...)
```

## Arguments

- atoms:

  a length-\\K\\ list where each element is an \\(N_k \times P)\\ matrix
  of atoms.

- marginals:

  marginal distributions for empirical measures; if `NULL` (default),
  uniform weights are set for all measures. Otherwise, it should be a
  length-\\K\\ list where each element is a length-\\N_i\\ vector of
  nonnegative weights that sum to 1.

- weights:

  weights for each individual measure; if `NULL` (default), each measure
  is considered equally. Otherwise, it should be a length-\\K\\ vector.

- num_support:

  the number of support points \\M\\ for the barycenter (default: 100).

- ...:

  extra parameters including

  abstol

  :   stopping criterion for iterations (default: 1e-6).

  maxiter

  :   maximum number of iterations (default: 10).

## Value

a list with three elements:

- support:

  an \\(M \times P)\\ matrix of the Wasserstein median's support points.

- weight:

  a length-\\M\\ vector of median's weights with all entries being
  \\1/M\\.

- history:

  a vector of cost values at each iteration.

## Examples

``` r
if (FALSE) { # \dontrun{
#-------------------------------------------------------------------
#     Free-Support Wasserstein Median of Multiple Gaussians
#
# * class 1 : samples from N((0,0),  Id)
# * class 2 : samples from N((20,0), Id)
#
#  We draw 8 empirical measures of size 50 from class 1, and 
#  2 from class 2. All measures have uniform weights.
#-------------------------------------------------------------------
## GENERATE DATA
#  8 empirical measures from class 1
input_measures = vector("list", length=10L)
for (i in 1:8){
  input_measures[[i]] = matrix(rnorm(50*2), ncol=2)
}
for (j in 9:10){
  base_draw = matrix(rnorm(50*2), ncol=2)
  base_draw[,1] = base_draw[,1] + 20
  input_measures[[j]] = base_draw
}

## COMPUTE
#  compute the Wasserstein median
run_median = rmedWB(input_measures, num_support = 50)
#  compute the Wasserstein barycenter
run_bary   = rbaryGD(input_measures, num_support = 50)

## VISUALIZE
opar <- par(no.readonly=TRUE)

#  draw the base points of two classes
base_1 = matrix(rnorm(80*2), ncol=2)
base_2 = matrix(rnorm(20*2), ncol=2)
base_2[,1] = base_2[,1] + 20
base_mat = rbind(base_1, base_2)
plot(base_mat, col="gray80", pch=19)

#  auxiliary information
title("estimated barycenter and median")
abline(v=0); abline(h=0)

#  draw the barycenter and the median
points(run_bary$support, col="red", pch=19)
points(run_median$support, col="blue", pch=19)
par(opar)
} # }
```
