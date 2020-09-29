#' Barycenter of Images by Cuturi & Doucet (2014)
#' 
#' Using entropic regularization for Wasserstein barycenter computation, \code{image14C} 
#' finds a \emph{barycentric} image \eqn{X^*} given multiple images \eqn{X_1,X_2,\ldots,X_N}. 
#' Please note the followings; (1) we only take a matrix as an image so please 
#' make it grayscale if not, (2) all images should be of same size - no resizing is performed. 
#' 
#' @param images a length-\eqn{N} list of same-size image matrices of size \eqn{(m\times n)}.
#' @param p an exponent for the order of the distance (default: 2).
#' @param weights a weight of each image; if \code{NULL} (default), uniform weight is set. Otherwise, 
#' it should be a length-\eqn{N} vector of nonnegative weights. 
#' @param lambda a regularization parameter; if \code{NULL} (default), a paper's suggestion 
#' would be taken, or it should be a nonnegative real number.
#' @param ... extra parameters including \describe{
#' \item{abstol}{stopping criterion for iterations (default: 1e-8).}
#' \item{init.image}{an initial weight image (default: uniform weight).}
#' \item{maxiter}{maximum number of iterations (default: 496).}
#' \item{nthread}{number of threads for OpenMP run (default: 1).}
#' \item{print.progress}{a logical to show current iteration (default: \code{TRUE}).}
#' }
#' 
#' @return an \eqn{(m\times n)} matrix of the barycentric image.
#' 
#' @references 
#' \insertRef{cuturi_fast_2014}{T4transport}
#' 
#' @examples 
#' #----------------------------------------------------------------------
#' #                       MNIST Data with Digit 3
#' #
#' # EXAMPLE 1 : Very Small  Example for CRAN; just showing how to use it!
#' # EXAMPLE 2 : Medium-size Example for Evolution of Output
#' #----------------------------------------------------------------------
#' # EXAMPLE 1
#' data(digit3)
#' datsmall = digit3[1:2]
#' outsmall = image14C(datsmall, maxiter=3)
#' 
#' \dontrun{
#' # EXAMPLE 2 : Barycenter of 100 Images
#' # RANDOMLY SELECT THE IMAGES
#' dat2 = digit3[sample(1:2000, 100)]  # select 100 images
#' 
#' # RUN SEQUENTIALLY
#' run10 = image14C(dat2, maxiter=10)                   # first 10 iterations
#' run20 = image14C(dat2, maxiter=10, init.image=run10) # run 40 more
#' run50 = image14C(dat2, maxiter=30, init.image=run20) # run 50 more
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3), pty="s")
#' image(dat2[[sample(100,1)]], axes=FALSE, main="a random image")
#' image(dat2[[sample(100,1)]], axes=FALSE, main="a random image")
#' image(dat2[[sample(100,1)]], axes=FALSE, main="a random image")
#' image(run10, axes=FALSE, main="barycenter after 10 iter")
#' image(run20, axes=FALSE, main="barycenter after 20 iter")
#' image(run50, axes=FALSE, main="barycenter after 50 iter")
#' par(opar)
#' }
#' 
#' @concept imagecenter
#' @export
image14C <- function(images, p=2, weights=NULL, lambda=NULL, ...){
  # CHECK THE INPUT
  name.f  = "image14C"
  check.f = check_images(images, name.f)
  
  # GRID AND TRANSFORM
  imgsize = dim(images[[1]])
  coordx  = seq(from=0, to=1, length.out=imgsize[2])
  coordy  = seq(from=1, to=0, length.out=imgsize[1])
  coords  = expand.grid(coordx, coordy)
  
  dxy    = as.matrix(stats::dist(coords)) 
  nimage = length(images)
  
  # OTHER INFORMATION
  myp = max(1, as.double(p))
  mymarginal = list()
  for (i in 1:nimage){
    mymarginal[[i]] = as.vector(t(images[[i]]))
  }
  myweights = valid_multiple_weight(weights, nimage, name.f)
  myweights = myweights/base::sum(myweights)
  if ((length(lambda)==0)&&(is.null(lambda))){
    mylambda = 60/(stats::median(dxy)^myp) # choice of the paper
  } else {
    mylambda = 1/max(100*.Machine$double.eps, as.double(lambda)) 
  }
  params = list(...)
  pnames = names(params)
  myiter = max(1, round(ifelse((("maxiter")%in%pnames), params$maxiter, 496)))
  mytol  = max(100*.Machine$double.eps, as.double(ifelse(("abstol"%in%pnames), params$abstol, 1e-8)))
  myshow = as.logical(ifelse(("print.progress"%in%pnames), params$print.progress, TRUE))
  mynthr = max(1, round(ifelse(("nthread"%in%pnames), params$nthread, 1))) # OpenMP Threads
  
  nsupport = base::nrow(coords)
  if ("init.image" %in% pnames){
    par_init = as.vector(t(params$init.image))
    par_init = par_init/base::sum(par_init)
    if ((length(par_init)!=nsupport)||(any(par_init < 0))){
      stop(paste0("* image14C : 'init.image' should be of matching size as other images with nonnegative values."))
    }
  } else {
    par_init = rep(1/nsupport, nsupport)
  }
  
  # RUN, WRAP, AND RETURN
  myoutput = image_barysinkhorn14(dxy, mymarginal, myweights, myp, mylambda, myiter, mytol, myshow, par_init, mynthr)
  myshaped = matrix(myoutput, imgsize[1], imgsize[2], byrow=TRUE)
  return(myshaped)
}

# library(Barycenter)
# data("three")
# image050 = three[1:50]
# image100 = three[1:100]
# image200 = three[1:200]
# 
# hey1 = imageC14(image050, maxiter=1000)
# hey2 = imageC14(image100, maxiter=10)
# hey3 = imageC14(image200, maxiter=10)

# par(mfrow=c(1,3))
# image(matrix(hey1, 28, 28, byrow=TRUE))
# image(matrix(hey2, 28, 28, byrow=TRUE))
# image(matrix(hey3, 28, 28, byrow=TRUE))


# x11()
# par(mfrow=c(5,5), pty="s")
# for (i in 1:25){
#   image(three[[i]])
# }
# 
# data("digit3")
# pdat  = digit3[1:10]
# out10 = image14C(pdat, maxiter=10, print.progress=TRUE)
# out20 = image14C(pdat, maxiter=10, print.progress=TRUE, init.image=out10)
# out50 = image14C(pdat, maxiter=30, print.progress=TRUE, init.image=out20)
# out100 = image14C(pdat, maxiter=50, print.progress=TRUE, init.image=out50) 
# out200 = image14C(pdat, maxiter=100, print.progress=TRUE, init.image=out100) 
# out500 = image14C(pdat, maxiter=300, print.progress=TRUE, init.image=out200)
# 
# par(mfrow=c(1,3), pty="s")
# image(out10,  main="after 10 iterations")
# image(out100, main="after 100 iterations")
# image(out500, main="after 500 iterations")
