# Barycenter of Images according to Benamou et al. (2015)

Using entropic regularization for Wasserstein barycenter computation,
`imagebary15B` finds a *barycentric* image \\X^\*\\ given multiple
images \\X_1,X_2,\ldots,X_N\\. Please note the followings; (1) we only
take a matrix as an image so please make it grayscale if not, (2) all
images should be of same size - no resizing is performed.

## Usage

``` r
imagebary15B(images, p = 2, weights = NULL, lambda = NULL, ...)
```

## Arguments

- images:

  a length-\\N\\ list of same-size image matrices of size \\(m\times
  n)\\.

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

  init.image

  :   an initial weight image (default: uniform weight).

  maxiter

  :   maximum number of iterations (default: 496).

  nthread

  :   number of threads for OpenMP run (default: 1).

  print.progress

  :   a logical to show current iteration (default: `TRUE`).

## Value

an \\(m\times n)\\ matrix of the barycentric image.

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
#                       MNIST Data with Digit 3
#
# EXAMPLE 1 : Very Small  Example for CRAN; just showing how to use it!
# EXAMPLE 2 : Medium-size Example for Evolution of Output
#----------------------------------------------------------------------
# EXAMPLE 1
data(digit3)
datsmall = digit3[1:2]
outsmall = imagebary15B(datsmall, maxiter=3)

if (FALSE) { # \dontrun{
# EXAMPLE 2 : Barycenter of 100 Images
# RANDOMLY SELECT THE IMAGES
data(digit3)
dat2 = digit3[sample(1:2000, 100)]  # select 100 images

# RUN SEQUENTIALLY
run05 = imagebary15B(dat2, maxiter=5)                    # first 5 iterations
run10 = imagebary15B(dat2, maxiter=5,  init.image=run05) # run 5 more
run50 = imagebary15B(dat2, maxiter=40, init.image=run10) # run 40 more

# VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,3), pty="s")
image(dat2[[sample(100,1)]], axes=FALSE, main="a random image")
image(dat2[[sample(100,1)]], axes=FALSE, main="a random image")
image(dat2[[sample(100,1)]], axes=FALSE, main="a random image")
image(run05, axes=FALSE, main="barycenter after 05 iter")
image(run10, axes=FALSE, main="barycenter after 10 iter")
image(run50, axes=FALSE, main="barycenter after 50 iter")
par(opar)
} # }
```
