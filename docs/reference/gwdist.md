# Gromov-Wasserstein Distance

Computes the Gromov-Wasserstein (GW) distance between two metric measure
spaces. Given two distance matrices \\D_x\\ and \\D_y\\ along with their
respective marginal distributions, the function solves the GW
optimization problem to obtain both the distance value and an associated
optimal transport plan.

The GW distance provides a way to compare datasets that may not lie in
the same ambient space by focusing on the intrinsic geometric structure
encoded in the pairwise distances. This implementation supports multiple
optimization schemes, including majorization–minimization (MM), proximal
gradient (PG), and Frank–Wolfe (FW).

## Usage

``` r
gwdist(Dx, Dy, wx = NULL, wy = NULL, ...)
```

## Arguments

- Dx:

  an \\(M\times M)\\ distance matrix or a
  [`dist`](https://rdrr.io/r/stats/dist.html) object of compatible size.

- Dy:

  an \\(N\times N)\\ distance matrix or a
  [`dist`](https://rdrr.io/r/stats/dist.html) object of compatible size.

- wx:

  a length-\\M\\ marginal density that sums to \\1\\. If `NULL`
  (default), uniform weight is set.

- wy:

  a length-\\N\\ marginal density that sums to \\1\\. If `NULL`
  (default), uniform weight is set.

- ...:

  extra parameters including

  maxiter

  :   maximum number of iterations (default: 496).

  abstol

  :   stopping criterion for iterations (default: 1e-10).

  method

  :   optimization method to use; can be one of `"mm"`, `"pg"`, or
      `"fw"` (default).

## Value

a named list containing

- distance:

  the computed GW distance value.

- plan:

  an \\(M\times N)\\ nonnegative matrix for the optimal transport plan.

## References

Mémoli F (2011). “Gromov–Wasserstein Distances and the Metric Approach
to Object Matching.” *Foundations of Computational Mathematics*,
**11**(4), 417–487. ISSN 1615-3375, 1615-3383.

## Examples

``` r
if (FALSE) { # \dontrun{
#-------------------------------------------------------------------
#                           Description
#
# * class 1 : iris dataset (columns 1-4) with perturbations
# * class 2 : class 1 rotated randomly in R^4
# * class 3 : samples from N((0,0), I)
#
#  We draw 10 empirical measures from each and compare 
#  the regular Wasserstein and GW distance. It is expected that 
#  the GW distance between class 1 and class 2 is negligible, 
#  while the regular Wasserstein distance is large. For simplicity, 
#  limit the cardinalities to 20.
#-------------------------------------------------------------------
## GENERATE DATA
set.seed(10)

#  prepare empty lists
inputs = vector("list", length=30)

#  generate class 1 and 2
iris_mat = as.matrix(iris[sample(1:150,20),1:4])
for (i in 1:10){
  inputs[[i]] = iris_mat + matrix(rnorm(20*4), ncol=4)
  inputs[[i+10]] = inputs[[i]]%*%qr.Q(qr(matrix(runif(16), ncol=4)))
}
#  generate class 3
for (j in 21:30){
  inputs[[j]] = matrix(rnorm(20*4), ncol=4)
}

## COMPUTE
#  empty arrays
dist_RW = array(0, c(30, 30))
dist_GW = array(0, c(30, 30))

#  compute pairwise distances
for (i in 1:29){
  X <- inputs[[i]]
  Dx <- stats::dist(X)
  for (j in (i+1):30){
  Y <- inputs[[j]]
  Dy <- stats::dist(Y)
  dist_RW[i,j] <- dist_RW[j,i] <- wasserstein(X, Y)$distance
  dist_GW[i,j] <- dist_GW[j,i] <- gwdist(Dx, Dy)$distance
  }
}

## VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2), pty="s")
image(dist_RW, xaxt="n", yaxt="n", main="Regular Wasserstein distance")
image(dist_GW, xaxt="n", yaxt="n", main="Gromov-Wasserstein distance")
par(opar)
} # }
```
