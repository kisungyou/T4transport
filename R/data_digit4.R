#' MNIST Images of Digit 4
#' 
#' \code{digit4} contains 2000 images from the famous MNIST dataset of digit 3. 
#' Each element of the list is an image represented as an \eqn{(28\times 28)} 
#' matrix that sums to 1. This normalization is conventional and it does not 
#' hurt its visualization via a basic `image()` function.
#' 
#' @usage data(digit4)
#' 
#' @format a length-\eqn{2000} named list \code{"digit4"} of \eqn{(28\times 28)} matrices.
#' 
#' @examples 
#' ## LOAD THE DATA
#' data(digit4)
#' 
#' ## SHOW A FEW
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,4), pty="s")
#' for (i in 1:8){
#'   image(digit4[[i]])
#' }
#' par(opar)
#' 
#' @concept data
"digit4"


# data is fetched by 'dslabs' package.
# id4 = which(mnist$train$labels==4)
# id4 = id4[sample(1:length(id4), 2000)]
# 
# digit4 = list()
# for (i in 1:2000){
#   digit4[[i]] = matrix(mnist$train$images[id4[i],], nrow=28)[,28:1]
# }
# for (i in 1:2000){
#   tgt = digit4[[i]]
#   digit4[[i]] = tgt/base::sum(as.vector(tgt))
# }
# 
# par(mfrow=c(5,5))
# for (i in 1:25){
#   image(digit4[[i]])
# }