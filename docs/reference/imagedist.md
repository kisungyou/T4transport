# Wasserstein Distance between Two Images

Given two grayscale images represented as numeric matrices, compute
their Wasserstein distance using an exact balanced optimal transport
solver. Each image is interpreted as a discrete probability distribution
on a common \\(m\times n)\\ grid. The ground cost is defined using the
Euclidean distance between grid locations.

## Usage

``` r
imagedist(x, y, p = 2)
```

## Arguments

- x:

  a grayscale image matrix of size \\(m\times n)\\ with nonnegative
  entries.

- y:

  a grayscale image matrix of size \\(m\times n)\\ with nonnegative
  entries.

- p:

  an exponent for the order of the distance (default: 2).

## Value

a list containing

- distance:

  the Wasserstein distance \\W_p(x,y)\\.

- plan:

  the optimal transport plan matrix of size \\(mn\times mn)\\.

## Examples

``` r
# \donttest{
#----------------------------------------------------------------------
#                       Small MNIST-like Example
#----------------------------------------------------------------------
# DATA
data(digit3)
x <- digit3[[1]]
y <- digit3[[2]]

# COMPUTE
W1 <- imagedist(x, y, p=1)
W2 <- imagedist(x, y, p=2)

# SHOW RESULTS
print(paste0("Wasserstein-1 distance: ", round(W1$distance,4)))
#> [1] "Wasserstein-1 distance: 0.0998"
print(paste0("Wasserstein-2 distance: ", round(W2$distance,4)))
#> [1] "Wasserstein-2 distance: 0.1267"
# }
```
