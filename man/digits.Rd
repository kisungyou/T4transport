% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_digits.R
\docType{data}
\name{digits}
\alias{digits}
\title{MNIST Images of All Digits}
\format{
a named list \code{"digits"} containing \describe{
\item{image}{length-5000 list of \eqn{(28\times 28)} image matrices.}
\item{label}{length-5000 vector of class labels from 0 to 9.}
}
}
\usage{
data(digits)
}
\description{
\code{digits} contains 5000 images from the famous MNIST dataset of all digits, 
consisting of 500 images per digit class from 0 to 9. 
Each digit image is represented as an \eqn{(28\times 28)} 
matrix that sums to 1. This normalization is conventional and it does not 
hurt its visualization via a basic `image()` function.
}
\examples{
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

}
\concept{data}
\keyword{datasets}
