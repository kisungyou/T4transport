#' Procrustes-Wasserstein Barycenter
#' 
#' @description
#' For a collection of empirical measures \eqn{\lbrace \mu_k\rbrace_{k=1}^K}, this 
#' function computes the Procrustes-Wasserstein (PW) barycenter \insertCite{adamo_2025_DepthLookProcrustesWasserstein}{T4transport}, 
#' which accounts for both measure transport and alignment 
#' through action of the orthogonal group.
#' 
#' @param atoms a length-\eqn{K} list where each element is an \eqn{(N_k \times P)} matrix of atoms.
#' @param marginals marginal distributions for empirical measures; if \code{NULL} (default), uniform weights are set for all measures. Otherwise, it should be a length-\eqn{K} list where each element is a length-\eqn{N_i} vector of nonnegative weights that sum to 1.
#' @param weights weights for each individual measure; if \code{NULL} (default), each measure is considered equally. Otherwise, it should be a length-\eqn{K} vector.
#' @param num_support the number of support points \eqn{M} for the PW barycenter (default: 100).
#' @param ... extra parameters including \describe{
#' \item{abstol}{stopping criterion for iterations (default: 1e-6).}
#' \item{maxiter}{maximum number of iterations (default: 10).}
#' }
#' 
#' @return a list with three elements:
#' \describe{
#'   \item{support}{an \eqn{(M \times P)} matrix of the PW barycenter's support points.}
#'   \item{weight}{a length-\eqn{M} vector of median's weights with all entries being \eqn{1/M}.}
#' }
#' 
#' @examples
#' \dontrun{
#' #-------------------------------------------------------------------
#' #         Free-Support PW Barycenter of Multiple Gaussians
#' #
#' # * class 1 : samples from N((0,0),  diag(c(4,1/4)))
#' # * class 2 : samples from N((10,0), diag(c(1/4,4)))
#' # * class 3 : samples from N((10,0), Id) randomly rotated
#' #
#' #  We draw 10 empirical measures from each and compare 
#' #  their barycenters under the regular and PW geometries.
#' #-------------------------------------------------------------------
#' ## GENERATE DATA
#' set.seed(10)
#' 
#' #  prepare empty lists
#' input_1 = vector("list", length=10L)
#' input_2 = vector("list", length=10L)
#' input_3 = vector("list", length=10L)
#' 
#' #  generate
#' random_rot = qr.Q(qr(matrix(runif(4), ncol=2)))
#' for (i in 1:10){
#'   input_1[[i]] = cbind(rnorm(50, sd=2), rnorm(50, sd=0.5))
#' }
#' for (j in 1:10){
#'   base_draw = cbind(rnorm(50, sd=0.5), rnorm(50, sd=2))
#'   base_draw[,1] = base_draw[,1] + 10
#'   
#'   input_2[[j]] = base_draw
#'   input_3[[j]] = base_draw%*%random_rot
#' }
#' 
#' ## COMPUTE
#' #  regular Wasserstein barycenters
#' regular_1 = rbaryGD(input_1, num_support=50)
#' regular_2 = rbaryGD(input_2, num_support=50)
#' regular_3 = rbaryGD(input_3, num_support=50)
#' 
#' #  Procrustes-Wasserstein barycenters
#' pw_1 = pwbary(input_1, num_support=50)
#' pw_2 = pwbary(input_2, num_support=50)
#' pw_3 = pwbary(input_3, num_support=50)
#' 
#' ## VISUALIZE
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(3,1))
#' 
#' #  set the x- and y-limits for display
#' lim_x = c(-12, 12)
#' lim_y = c(-10, 5)
#' 
#' #  plot prototypical measures per class
#' plot(input_1[[1]], pch=19, cex=0.5, col="gray80", 
#'      main="3 types of measures", xlab="", ylab="",
#'      xlim=lim_x, ylim=lim_y)
#' points(input_2[[1]], pch=19, cex=0.5, col="gray50")
#' points(input_3[[1]], pch=19, cex=0.5, col="gray10")
#' 
#' #  plot regular barycenters
#' plot(regular_1$support, pch=19, cex=0.5, col="blue", 
#'      main="Regular Wasserstein barycenters",
#'      xlab="", ylab="", xlim=lim_x, ylim=lim_y)
#' points(regular_2$support, pch=19, cex=0.5, col="cyan")
#' points(regular_3$support, pch=19, cex=0.5, col="red")
#' 
#' #  plot PW barycenters
#' plot(pw_1$support, pch=19, cex=0.5, col="blue", 
#'      main="Procrustes-Wasserstein barycenters",
#'      xlab="", ylab="", xlim=lim_x, ylim=lim_y)
#' points(pw_2$support, pch=19, cex=0.5, col="cyan")
#' points(pw_3$support, pch=19, cex=0.5, col="red")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{adamo_2025_DepthLookProcrustesWasserstein}{T4transport}
#' 
#' @concept free_centroid
#' @export
pwbary <- function(atoms, marginals=NULL, weights=NULL, num_support=100, ...){
  ## INPUT : EXPLICIT
  name.f    = "pwbary"
  par_measures   = valid_multiple_measures(atoms, base::ncol(atoms[[1]]), name.f)
  num_atoms      = unlist(lapply(atoms, nrow))
  par_marginals  = valid_multiple_marginal(marginals, num_atoms, name.f)
  par_weights    = valid_multiple_weight(weights, length(atoms), name.f)
  par_numsupport = max(2, round(num_support))
  
  ## INPUT : IMPLICIT
  params = list(...)
  pnames = names(params)
  
  if ("maxiter" %in% pnames){
    par_maxiter = max(1, round(params[["maxiter"]]))
  } else {
    par_maxiter = 10
  }
  if ("abstol" %in% pnames){
    par_abstol = max(10*.Machine$double.eps, as.double(params[["abstol"]]))
  } else {
    par_abstol = 1e-6
  }
  
  ## COMPUTE
  #  initialize the measure
  init_measure = aux_ginit(par_measures, par_numsupport)
  
  #  compute with the internal routine
  support_out = pwbary_internal(
    X_init = init_measure,
    Xs = par_measures,
    ws = par_marginals,
    rs = par_weights,
    maxiter = par_maxiter,
    abstol = par_abstol
  )
  
  # return
  output = list()
  output[["support"]] = support_out
  output[["weight"]] = rep(1/par_numsupport, par_numsupport)
  return(output)
}




# internal routine for pwbary ---------------------------------------------
# X_init: initial support matrix
# Xs: list of input supports
# ws: list of input weights
# rs: vector of relative weights
# maxiter: maximum number of iterations
# abstol: absolute tolerance
#' @keywords internal
#' @noRd
pwbary_internal <- function(X_init, Xs, ws, rs, maxiter, abstol){
  # initialize
  N = length(Xs)
  X_old = X_init
  X_num = base::nrow(X_old)
  X_weight = rep(1/X_num, X_num)
  
  # iterate
  for (it in seq_len(round(maxiter))){
    # step 1 : compute optimal plan and alignment
    opt_plans = vector("list", length=N)
    opt_align = vector("list", length=N)
    for (n in 1:N){
      compute_now = cpp_pwdist(X_old, X_weight, Xs[[n]], ws[[n]], maxiter, abstol)
      opt_plans[[n]] = as.matrix(compute_now[["est_plan"]])
      opt_align[[n]] = as.matrix(compute_now[["est_P"]])
    }
    
    # step 2 : update the support
    X_new = array(0, c(X_num, ncol(X_old)))
    for (n in 1:N){
      X_new = X_new + rs[n]*(diag(1/X_weight)%*%(opt_plans[[n]]%*%(Xs[[n]]%*%opt_align[[n]])))
    }
    
    # stop
    X_inc = base::norm(X_new-X_old,"F")
    X_old = X_new
    if (X_inc < abstol){
      break
    }
  }
  
  # return
  return(X_old)
}
# apply(rbind(apply(regular_1$support, 2, min),
#       apply(regular_2$support, 2, min),
#       apply(regular_3$support, 2, min)), 2, min)
# 
# apply(rbind(apply(regular_1$support, 2, max),
#       apply(regular_2$support, 2, max),
#       apply(regular_3$support, 2, max)), 2, max)