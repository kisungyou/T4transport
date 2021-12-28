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
  imgsize = dim(images[[1]])
  coordx  = seq(from=0, to=1, length.out=imgsize[2])
  coordy  = seq(from=1, to=0, length.out=imgsize[1])
  coords  = expand.grid(coordx, coordy)
  
  dxy    = as.matrix(stats::dist(coords)) 
  nimage = length(images)
  
  # --------------------------------------------------------------------------
  # OTHER INFORMATION
  myp = 2
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
  mynthr = max(1, round(ifelse(("nthread"%in%pnames), params$nthread, 1))) # OpenMP Threads
  

  
  nsupport = base::nrow(coords)
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
  digit_show = abs(round(log10(mytol)))+1
  
  # --------------------------------------------------------------------------
  # ITERATIVE PROCEDURE
  # RUN, WRAP, AND RETURN
  image_old = par_init
  for (it in 1:myiter){
    # update the relative weight
    for (i in 1:length(myweights)){
      sinkhorn_run = cpp_sinkhorn13(as.vector(mymarginal[[i]]), image_old, dxy, mylambda, myp, myiter, mytol)
      myweights[i] = 1/as.double(sinkhorn_run$distance) # compute distance by Sinkhorn
    }
    myweights = myweights/base::sum(myweights)
    
    # update an image
    image_new = routine_bary15B(dxy, mymarginal, myweights, myp, mylambda, myiter, mytol, FALSE, image_old, mynthr)
    image_new = as.vector(image_new)
    
    # compute the error & update
    increment = sqrt(sum(image_old-image_new)^2)
    if (increment < mytol){
      if (myshow){
        print(paste0("* medimage22Y : algorithm terminates at iteration ",it," : increment error=",round(increment,digit_show)))
      }
      return(matrix(image_new, imgsize[1], imgsize[2], byrow=TRUE))
      break
    }
    image_old = image_new
    if (myshow){
      print(paste0("* medimage22Y : iteration ",it,"/",myiter," complete : increment error=",round(increment,digit_show)))
    }
  }
  return(matrix(image_old, imgsize[1], imgsize[2], byrow=TRUE))
}