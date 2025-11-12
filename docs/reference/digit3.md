# MNIST Images of Digit 3

`digit3` contains 2000 images from the famous MNIST dataset of digit 3.
Each element of the list is an image represented as an \\(28\times 28)\\
matrix that sums to 1. This normalization is conventional and it does
not hurt its visualization via a basic \`image()\` function.

## Usage

``` r
data(digit3)
```

## Format

a length-\\2000\\ named list `"digit3"` of \\(28\times 28)\\ matrices.

## Examples

``` r
## LOAD THE DATA
data(digit3)

## SHOW A FEW
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,4), pty="s")
for (i in 1:8){
  image(digit3[[i]])
}

par(opar)
```
