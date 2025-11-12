# Barycenter of Histograms by Benamou et al. (2015)

Given multiple histograms represented as `"histogram"` S3 objects,
compute Wasserstein barycenter. We need one requirement that all
histograms in an input list `hists` must have **same breaks**. See the
example on how to construct a histogram on predefined breaks/bins.

## Usage

``` r
histbary15B(hists, p = 2, weights = NULL, lambda = NULL, ...)
```

## Arguments

- hists:

  a length-\\N\\ list of histograms (`"histogram"` object) of same
  breaks.

- p:

  an exponent for the order of the distance (default: 2).

- weights:

  a weight of each image; if `NULL` (default), uniform weight is set.
  Otherwise, it should be a length-\\N\\ vector of nonnegative weights.

- lambda:

  a regularization parameter; if `NULL` (default), a paper's suggestion
  would be taken, or it should be a nonnegative real number.

- ...:

  extra parameters including

  abstol

  :   stopping criterion for iterations (default: 1e-8).

  init.vec

  :   an initial weight vector (default: uniform weight).

  maxiter

  :   maximum number of iterations (default: 496).

  nthread

  :   number of threads for OpenMP run (default: 1).

  print.progress

  :   a logical to show current iteration (default: `TRUE`).

## Value

a `"histogram"` object of barycenter.

## References

Benamou J, Carlier G, Cuturi M, Nenna L, Peyré G (2015). “Iterative
Bregman Projections for Regularized Transportation Problems.” *SIAM
Journal on Scientific Computing*, **37**(2), A1111-A1138. ISSN
1064-8275, 1095-7197,
[doi:10.1137/141000439](https://doi.org/10.1137/141000439) .

## See also

[`fbary15B`](https://www.kisungyou.com/T4transport/reference/fbary15B.md)

## Examples

``` r
#----------------------------------------------------------------------
#                      Binned from Two Gaussians
#
# EXAMPLE : Very Small Example for CRAN; just showing how to use it!
#----------------------------------------------------------------------
# GENERATE FROM TWO GAUSSIANS WITH DIFFERENT MEANS
set.seed(100)
x  = stats::rnorm(1000, mean=-4, sd=0.5)
y  = stats::rnorm(1000, mean=+4, sd=0.5)
bk = seq(from=-10, to=10, length.out=20)

# HISTOGRAMS WITH COMMON BREAKS
histxy = list()
histxy[[1]] = hist(x, breaks=bk, plot=FALSE)
histxy[[2]] = hist(y, breaks=bk, plot=FALSE)

# COMPUTE
hh = histbary15B(histxy, maxiter=5)

# VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(1,2))
barplot(histxy[[1]]$density, col=rgb(0,0,1,1/4), 
        ylim=c(0, 0.75), main="Two Histograms")
barplot(histxy[[2]]$density, col=rgb(1,0,0,1/4), 
        ylim=c(0, 0.75), add=TRUE)
barplot(hh$density, main="Barycenter",
        ylim=c(0, 0.75))

par(opar)
```
