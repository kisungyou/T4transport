#' Barycenter of Images according to Benamou et al. (2015)
#' 
#' Using entropic regularization for Wasserstein barycenter computation, \code{imagebary15B} 
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
#' outsmall = imagebary15B(datsmall, maxiter=3)
#' 
#' \dontrun{
#' # EXAMPLE 2 : Barycenter of 100 Images
#' # RANDOMLY SELECT THE IMAGES
#' data(digit3)
#' dat2 = digit3[sample(1:2000, 100)]  # select 100 images
#' 
#' # RUN SEQUENTIALLY
#' run05 = imagebary15B(dat2, maxiter=5)                    # first 5 iterations
#' run10 = imagebary15B(dat2, maxiter=5,  init.image=run05) # run 5 more
#' run50 = imagebary15B(dat2, maxiter=40, init.image=run10) # run 40 more
#' 
#' # VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3), pty="s")
#' image(dat2[[sample(100,1)]], axes=FALSE, main="a random image")
#' image(dat2[[sample(100,1)]], axes=FALSE, main="a random image")
#' image(dat2[[sample(100,1)]], axes=FALSE, main="a random image")
#' image(run05, axes=FALSE, main="barycenter after 05 iter")
#' image(run10, axes=FALSE, main="barycenter after 10 iter")
#' image(run50, axes=FALSE, main="barycenter after 50 iter")
#' par(opar)
#' }
#' 
#' @seealso \code{\link{fbary15B}}
#' 
#' @references 
#' \insertRef{benamou_2015_IterativeBregmanProjections}{T4transport}
#' 
#' @concept image
#' @export
imagebary15B <- function(images, p=2, weights=NULL, lambda=NULL, ...){
  # CHECK THE INPUT
  name.f  = "imagebary15B"
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
    mylambda = 1/(60/(stats::median(dxy)^myp)) # choice of the paper
  } else {
    mylambda = max(100*.Machine$double.eps, as.double(lambda)) 
  }
  params = list(...)
  pnames = names(params)
  myiter = max(1, round(ifelse((("maxiter")%in%pnames), params$maxiter, 496)))
  mytol  = max(100*.Machine$double.eps, as.double(ifelse(("abstol"%in%pnames), params$abstol, 1e-8)))
  myshow = as.logical(ifelse(("print.progress"%in%pnames), params$print.progress, FALSE))
  mynthr = max(1, round(ifelse(("nthread"%in%pnames), params$nthread, 1))) # OpenMP Threads
  
  nsupport = base::nrow(coords)
  if ("init.image" %in% pnames){
    par_init = as.vector(t(params$init.image))
    par_init = par_init/base::sum(par_init)
    if ((length(par_init)!=nsupport)||(any(par_init < 0))){
      stop(paste0("* imagebary15B : 'init.image' should be of matching size as other images with nonnegative values."))
    }
  } else {
    par_init = rep(1/nsupport, nsupport)
  }
  
  # RUN, WRAP, AND RETURN
  myoutput = routine_bary15B(dxy, mymarginal, myweights, myp, mylambda, myiter, mytol, myshow, par_init, mynthr)
  myshaped = matrix(myoutput, imgsize[1], imgsize[2], byrow=TRUE)
  return(myshaped)
}
