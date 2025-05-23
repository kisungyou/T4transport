% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/image_bary15B.R
\name{imagebary15B}
\alias{imagebary15B}
\title{Barycenter of Images according to Benamou et al. (2015)}
\usage{
imagebary15B(images, p = 2, weights = NULL, lambda = NULL, ...)
}
\arguments{
\item{images}{a length-\eqn{N} list of same-size image matrices of size \eqn{(m\times n)}.}

\item{p}{an exponent for the order of the distance (default: 2).}

\item{weights}{a weight of each image; if \code{NULL} (default), uniform weight is set. Otherwise, 
it should be a length-\eqn{N} vector of nonnegative weights.}

\item{lambda}{a regularization parameter; if \code{NULL} (default), a paper's suggestion 
would be taken, or it should be a nonnegative real number.}

\item{...}{extra parameters including \describe{
\item{abstol}{stopping criterion for iterations (default: 1e-8).}
\item{init.image}{an initial weight image (default: uniform weight).}
\item{maxiter}{maximum number of iterations (default: 496).}
\item{nthread}{number of threads for OpenMP run (default: 1).}
\item{print.progress}{a logical to show current iteration (default: \code{TRUE}).}
}}
}
\value{
an \eqn{(m\times n)} matrix of the barycentric image.
}
\description{
Using entropic regularization for Wasserstein barycenter computation, \code{imagebary15B} 
finds a \emph{barycentric} image \eqn{X^*} given multiple images \eqn{X_1,X_2,\ldots,X_N}. 
Please note the followings; (1) we only take a matrix as an image so please 
make it grayscale if not, (2) all images should be of same size - no resizing is performed.
}
\examples{
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

\dontrun{
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
}

}
\references{
\insertRef{benamou_2015_IterativeBregmanProjections}{T4transport}
}
\seealso{
\code{\link{fbary15B}}
}
\concept{image}
