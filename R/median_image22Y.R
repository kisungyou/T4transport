#' Wasserstein Median of Images by You et al. (2022)
#' 
#' 
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
#' datsmall = digit3[1:3]
#' outsmall = medimage22Y(datsmall, print.progress=TRUE)
#' 
#' 
#' @concept median_image
#' @export
medimage22Y <- function(images, weights=NULL, lambda=NULL, ...){
  # --------------------------------------------------------------------------
  # CHECK THE INPUT
  name.f  = "medimage22Y"
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
      stop(paste0("* medimage22Y : 'init.image' should be of matching size as other images with nonnegative values."))
    }
  } else {
    par_init = rep(1/nsupport, nsupport)
  }
  if ("print.progress"%in%pnames){
    myshow = as.logical(params$print.progress)
  } else {
    myshow = FALSE
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
        print(paste0("* medimage22Y : algorithm terminates at iteration ",it," : increment error=",increment))
      }
      return(matrix(image_new, imgsize[1], imgsize[2], byrow=TRUE))
      break
    }
    image_old = image_new
    if (myshow){
      print(paste0("* medimage22Y : iteration ",it,"/",myiter," complete : increment error=",increment))
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