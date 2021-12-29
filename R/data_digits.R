#' MNIST Images of All Digits
#' 
#' \code{digits} contains 5000 images from the famous MNIST dataset of all digits, 
#' consisting of 500 images per digit class from 0 to 9. 
#' Each digit image is represented as an \eqn{(28\times 28)} 
#' matrix that sums to 1. This normalization is conventional and it does not 
#' hurt its visualization via a basic `image()` function.
#' 
#' @usage data(digits)
#' 
#' @format a named list \code{"digits"} containing \describe{
#' \item{image}{length-5000 list of \eqn{(28\times 28)} image matrices.}
#' \item{label}{length-5000 vector of class labels from 0 to 9.}
#' }
#' 
#' @examples 
#' ## LOAD THE DATA
#' data(digits)
#' 
#' ## SHOW A FEW
#' #  Select 9 random images
#' subimgs = digits$image[sample(1:5000, 9)]
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(3,3), pty="s")
#' for (i in 1:9){
#'   image(subimgs[[i]])
#' }
#' par(opar)
#' 
#' @concept data
"digits"

# # data is fetched by 'dslabs' package.
# library(dslabs)
# mnist=dslabs::read_mnist()
# 
# ulabel  = sort(unique(mnist$train$labels))
# labels  = c()
# images  = c()
# counter = 0
# for (i in 0:9){
#   id_now  = which(mnist$train$labels==i)
#   im_now  = id_now[sample(1:length(id_now),500)]
#   
#   labels = c(labels, rep(i, 500))
#   for (j in 1:500){
#     counter = counter + 1
#     tgt     = matrix(mnist$train$images[im_now[j],], nrow=28)[,28:1]
#     images[[counter]] = tgt/base::sum(tgt)
#   }
# }
# 
# digits = list(image=images, label=labels)
# usethis::use_data(digits)