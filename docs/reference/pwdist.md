# Procrustes-Wasserstein Distance

Given two empirical measures \$\$\mu = \sum\_{m=1}^M \mu_m
\delta\_{X_m}\quad\textrm{and}\quad \nu = \sum\_{n=1}^N \nu_n
\delta\_{Y_n}\$\$ in \\\mathbb{R}^P\\, the Procrustes-Wasserstein (PW)
distance is defined as follows: \$\$ PW_2^2(\mu, \nu) = \min\_{Q\in
\mathcal{O}(P)} W_2^2(\mu, Q\_\\ \nu), \$\$ where \\\mathcal{O}(P)\\ is
the orthogonal group and \\Q\_\\\nu\\ is the pushforward via \\Q\\.

## Usage

``` r
pwdist(X, Y, wx = NULL, wy = NULL, ...)
```

## Arguments

- X:

  an \\(M\times P)\\ matrix of row observations.

- Y:

  an \\(N\times P)\\ matrix of row observations.

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

## Value

a named list containing

- distance:

  the computed PW distance value.

- plan:

  an \\(M\times N)\\ nonnegative matrix for the optimal transport plan.

- alignment:

  an optimal alignment matrix of size \\(P\times P)\\ in
  \\\mathcal{O}(P)\\.

## References

Adamo D, Corneli M, Vuillien M, Vila E (2025). “An in Depth Look at the
Procrustes-Wasserstein Distance: Properties and Barycenters.” In
*Forty-Second International Conference on Machine Learning*.

## Examples

``` r
if (FALSE) { # \dontrun{
#-------------------------------------------------------------------
#                           Description
#
# * class 1 : samples from N((0,0),  diag(c(4,1/4)))
# * class 2 : samples from N((10,0), diag(c(1/4,4)))
# * class 3 : samples from N((10,0), diag(c(1/4,4))) randomly rotated
#
#  We draw 10 empirical measures from each and compare 
#  the regular Wasserstein and PW distance.
#-------------------------------------------------------------------
## GENERATE DATA
set.seed(10)

#  prepare empty lists
inputs = vector("list", length=30)

#  generate
random_rot = qr.Q(qr(matrix(runif(4), ncol=2)))
for (i in 1:10){
  inputs[[i]] = matrix(rnorm(50*2), ncol=2)
}
for (j in 11:20){
  base_draw = matrix(rnorm(50*2), ncol=2)
  base_draw[,1] = base_draw[,1] + 10
  
  inputs[[j]] = base_draw
  inputs[[j+10]] = base_draw%*%random_rot
}

## COMPUTE
#  empty arrays
dist_RW = array(0, c(30, 30))
dist_PW = array(0, c(30, 30))

#  compute pairwise distances
for (i in 1:29){
  for (j in (i+1):30){
  dist_RW[i,j] <- dist_RW[j,i] <- wasserstein(inputs[[i]], inputs[[j]])$distance
  dist_PW[i,j] <- dist_PW[j,i] <- pwdist(inputs[[i]], inputs[[j]])$distance
  }
}

## VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2), pty="s")
image(dist_RW, xaxt="n", yaxt="n", main="Regular Wasserstein distance")
image(dist_PW, xaxt="n", yaxt="n", main="PW distance")
par(opar)
} # }
```
