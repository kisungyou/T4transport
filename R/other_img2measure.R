#' Extract a discrete measure from a gray-scale image matrix
#' 
#' This function takes a gray-scale image represented as a matrix \eqn{X} and
#' converts it into a discrete measure suitable for optimal transport computations
#' in a Lagrangian framework. Pixel intensities are normalized to sum to one, and
#' the nonzero pixels are represented as weighted points (support and weights). 
#' 
#' @param X An \eqn{(N,2)} nonnegative matrix representing a gray-scale image, where each entry 
#' corresponds to a pixel intensity.
#' @param threshold A logical flag indicating whether to threshold very small weights smaller than machine epsilon.
#' 
#' @return A named list containing\describe{
#' \item{support}{an \eqn{(M\times 2)} matrix of coordinates for the nonzero pixels, where each row is a point \eqn{(x,y)}.}
#' \item{weight}{a length-\eqn{M} vector of weights corresponding to the nonzero pixels, summing to \eqn{1}.}
#' }
#' 
#' @examples
#' \donttest{
#' #-------------------------------------------------------------------
#' #                           Description
#' #
#' # Take a digit image and compare visualization.
#' #-------------------------------------------------------------------
#' # load the data and select the first image
#' data(digit3)
#' img_matrix = digit3[[1]]
#' 
#' # extract a discrete measure
#' img_measure = img2measure(img_matrix, threshold=TRUE)
#' w  <- img_measure$weight
#' w_norm <- w / max(w)          # now runs from 0 to 1
#' col_scale <- gray(1 - w_norm) # 1 = white, 0 = black
#' 
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(img_matrix, xaxt="n", yaxt="n", main="Image Matrix")
#' plot(img_measure$support, 
#'      col = col_scale, xlab="", ylab="",
#'      pch = 19, cex = 0.5, xaxt = "n", yaxt = "n",
#'      main = "Extracted Discrete Measure")
#' par(opar)
#' }
#' 
#' @concept other
#' @export
img2measure <- function(X, threshold=TRUE){
  # check
  if (!is.matrix(X)){
    stop("* img2measure: input 'X' should be a matrix.")
  }
  if (any(is.na(X))||any(is.infinite(X))||any(X < 0)){
    stop("* img2measure: input 'X' should have nonnegative finite entries only.")
  }
  
  # normalize into a measure
  Xnorm = X/base::sum(X)
  if (threshold){
    Xnorm[Xnorm < .Machine$double.eps] = 0
    Xnorm = Xnorm/base::sum(Xnorm) 
  }
  
  nr <- nrow(Xnorm)
  nc <- ncol(Xnorm)
  
  # indices of nonzero pixels
  nz_idx <- which(Xnorm > 0, arr.ind = TRUE)
  if (nrow(nz_idx) == 0L) {
    return(list(
      support = matrix(numeric(0), ncol = 2,
                       dimnames = list(NULL, c("x", "y"))),
      weight  = numeric(0)
    ))
  }
  
  w <- Xnorm[nz_idx]
  
  # --- coordinate system to match image() ---
  # x = row index (left→right), y = column index (bottom→top)
  # use pixel centers: index - 0.5
  support <- cbind(
    x = nz_idx[, "row"] - 0.5,
    y = nz_idx[, "col"] - 0.5
  )
  
  colnames(support) <- c("x", "y")
  
  return(list(
    support = support,
    weight  = w
  ))
}