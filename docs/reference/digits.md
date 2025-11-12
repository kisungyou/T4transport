# MNIST Images of All Digits

`digits` contains 5000 images from the famous MNIST dataset of all
digits, consisting of 500 images per digit class from 0 to 9. Each digit
image is represented as an \\(28\times 28)\\ matrix that sums to 1. This
normalization is conventional and it does not hurt its visualization via
a basic \`image()\` function.

## Usage

``` r
data(digits)
```

## Format

a named list `"digits"` containing

- image:

  length-5000 list of \\(28\times 28)\\ image matrices.

- label:

  length-5000 vector of class labels from 0 to 9.

## Examples

``` r
## LOAD THE DATA
data(digits)

## SHOW A FEW
#  Select 9 random images
subimgs = digits$image[sample(1:5000, 9)]

opar <- par(no.readonly=TRUE)
par(mfrow=c(3,3), pty="s")
for (i in 1:9){
  image(subimgs[[i]])
}

par(opar)
```
