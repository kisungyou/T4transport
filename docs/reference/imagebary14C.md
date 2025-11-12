# Barycenter of Images according to Cuturi & Doucet (2014)

Using entropic regularization for Wasserstein barycenter computation,
`imagebary14C` finds a *barycentric* image \\X^\*\\ given multiple
images \\X_1,X_2,\ldots,X_N\\. Please note the followings; (1) we only
take a matrix as an image so please make it grayscale if not, (2) all
images should be of same size - no resizing is performed.

## Usage

``` r
imagebary14C(images, p = 2, weights = NULL, lambda = NULL, ...)
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

Cuturi M, Doucet A (2014-06-22/2014-06-24). “Fast Computation of
Wasserstein Barycenters.” In Xing EP, Jebara T (eds.), *Proceedings of
the 31st International Conference on Machine Learning*, volume 32 of
*Proceedings of Machine Learning Research*, 685–693.

## See also

[`fbary14C`](https://www.kisungyou.com/T4transport/reference/fbary14C.md)

## Examples

``` r
if (FALSE) { # \dontrun{
#----------------------------------------------------------------------
#                       MNIST Data with Digit 3
#
# EXAMPLE 1 : Very Small  Example for CRAN; just showing how to use it!
# EXAMPLE 2 : Medium-size Example for Evolution of Output
#----------------------------------------------------------------------
# EXAMPLE 1
data(digit3)
datsmall = digit3[1:2]
outsmall = imagebary14C(datsmall, maxiter=3)

# EXAMPLE 2 : Barycenter of 100 Images
# RANDOMLY SELECT THE IMAGES
data(digit3)
dat2 = digit3[sample(1:2000, 100)]  # select 100 images

# RUN SEQUENTIALLY
run10 = imagebary14C(dat2, maxiter=10)                   # first 10 iterations
run20 = imagebary14C(dat2, maxiter=10, init.image=run10) # run 40 more
run50 = imagebary14C(dat2, maxiter=30, init.image=run20) # run 50 more

# VISUALIZE
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,3), pty="s")
image(dat2[[sample(100,1)]], axes=FALSE, main="a random image")
image(dat2[[sample(100,1)]], axes=FALSE, main="a random image")
image(dat2[[sample(100,1)]], axes=FALSE, main="a random image")
image(run10, axes=FALSE, main="barycenter after 10 iter")
image(run20, axes=FALSE, main="barycenter after 20 iter")
image(run50, axes=FALSE, main="barycenter after 50 iter")
par(opar)
} # }
```
