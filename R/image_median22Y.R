#' Wasserstein Median of Images by You et al. (2022)
#' 
#' Given multiple images \eqn{X_1,\ldots,X_N}, the Wasserstein median of 
#' order 2 is computed. The proposed method relies on a choice of barycenter computation 
#' in that we opt for an algorithm of \code{\link{imagebary15B}}, which uses 
#' entropic regularization for barycenter computation. Please note the followings; (1) we only take a matrix as an image so please 
#' make it grayscale if not, (2) all images should be of same size - no resizing is performed. 
#' 
#' @param images a length-\eqn{N} list of same-size image matrices of size \eqn{(m\times n)}.
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
#' @return an \eqn{(m\times n)} matrix of the Wasserstein median image.
#' 
#' @examples
#' \dontrun{
#' #----------------------------------------------------------------------
#' #                       MNIST Data with Digit 3
#' #
#' # EXAMPLE : Very Small Example for CRAN; just showing how to use it!
#' #----------------------------------------------------------------------
#' # LOAD THE DATA
#' data(digit3)
#' datsmall = digit3[1:10]
#'  
#' # COMPUTE
#' outsmall = imagemed22Y(datsmall, maxiter=5)
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,4), pty="s")
#' image(outsmall, xaxt='n', yaxt='n', main="Wasserstein Median")
#' image(datsmall[[3]], xaxt='n', yaxt='n', main="3rd image")
#' image(datsmall[[6]], xaxt='n', yaxt='n', main="6th image")
#' image(datsmall[[9]], xaxt='n', yaxt='n', main="9th image")
#' par(opar)
#' } 
#' 
#' @concept image
#' @export
imagemed22Y <- function(images, weights=NULL, lambda=NULL, ...){
  # --------------------------------------------------------------------------
  # CHECK THE INPUT
  name.f  = "imagemed22Y"
  check.f = check_images(images, name.f)
  
  # --------------------------------------------------------------------------
  # GRID AND TRANSFORM
  imgsize  = dim(images[[1]])
  coordx   = seq(from=0, to=1, length.out=imgsize[2])
  coordy   = seq(from=1, to=0, length.out=imgsize[1])
  coords   = expand.grid(coordx, coordy)
  nsupport = base::nrow(coords)
  
  dxy    = as.matrix(stats::dist(coords)) 
  nimage = length(images)
  
  # --------------------------------------------------------------------------
  # OTHER INFORMATION
  myp = 2
  mymarginal = list()
  for (i in 1:nimage){
    mymarginal[[i]] = as.vector(t(images[[i]]))
  }
  mypi      = valid_multiple_weight(weights, nimage, name.f)
  mypi      = mypi/base::sum(mypi)
  myweights = mypi
  if ((length(lambda)==0)&&(is.null(lambda))){
    mylambda = 1/(60/(stats::median(dxy)^myp)) # choice of the paper
  } else {
    mylambda = max(100*.Machine$double.eps, as.double(lambda)) 
  }
  
  params = list(...)
  pnames = names(params)
  
  if ("init.image" %in% pnames){
    par_init = as.vector(t(params$init.image))
    par_init = par_init/base::sum(par_init)
    if ((length(par_init)!=nsupport)||(any(par_init < 0))){
      stop(paste0("* imagemed22Y : 'init.image' should be of matching size as other images with nonnegative values."))
    }
  } else {
    par_init = rep(1/nsupport, nsupport)
  }
  if ("print.progress"%in%pnames){
    myshow = as.logical(params$print.progress)
  } else {
    myshow = TRUE
  }
  if ("maxiter"%in%pnames){
    myiter = max(1, round(params$maxiter))
  } else {
    myiter = 496
  }
  if ("abstol"%in%pnames){
    mytol = max(100*.Machine$double.eps, params$abstol)
  } else {
    mytol = 1e-8
  }
  round_digit = ceiling(abs(log10(mytol)))
  
  if ("nthread"%in%pnames){ # OpenMP Threads
    mynthr = max(1, round(params$nthread))
  } else {
    mynthr = 1
  }
  
  # --------------------------------------------------------------------------
  # ITERATIVE PROCEDURE
  # RUN, WRAP, AND RETURN
  image_old = par_init
  for (it in 1:myiter){
    # update the relative weight
    for (i in 1:length(myweights)){
      sinkhorn_run = cpp_sinkhorn13(as.vector(mymarginal[[i]]), image_old, dxy, mylambda, myp, 100, mytol)
      myweights[i] = mypi[i]/as.double(sinkhorn_run$distance) # compute distance by Sinkhorn
    }
    myweights = myweights/base::sum(myweights)
    
    # update an image
    image_new = routine_bary15B(dxy, mymarginal, myweights, myp, mylambda, 100, mytol, FALSE, image_old, mynthr)
    image_new = as.vector(image_new)
    
    # compute the error & update
    increment = max(as.vector(abs(image_old-image_new)))
    if (increment < mytol){
      if (myshow){
        print(paste0("* imagemed22Y : algorithm terminates at iteration ",it," : increment=",round(increment, round_digit)))
      }
      return(matrix(image_new, imgsize[1], imgsize[2], byrow=TRUE))
      break
    }
    image_old = image_new
    if (myshow){
      print(paste0("* imagemed22Y : iteration ",it,"/",myiter," complete : increment=",round(increment, round_digit)))
    }
  }
  return(matrix(image_old, imgsize[1], imgsize[2], byrow=TRUE))
}

# data("digit3")
# data("digit4")
# 
# gathered = c(digit3[sample(1:2000, 7)], digit4[sample(1:2000, 3)])
# fbary1   = image15B(gathered, p=1, maxiter=50, print.progress=TRUE)
# fbary2   = image15B(gathered, p=2, maxiter=50, print.progress=TRUE)
# fmedian  = medimage22Y(gathered, print.progress=TRUE, maxiter=20)
# 
# par(mfrow=c(1,3), pty="s")
# image(fbary1)
# image(fbary2)
# image(fmedian)
# 
# par(mfrow=c(2,5), pty="s")
# for (i in 1:10){
#   image(gathered[[i]])
# }
