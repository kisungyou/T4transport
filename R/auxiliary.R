## AUXILIARY FUNCTIONS
#  (1) valid_weight            : weight vector of provided length
#  (2) valid_distance          : distance matrix should contain no {(-)/Inf,NA}s
#  (3) valid_matrixed          : vector -> matrix.
#  (4) valid_multiple_measures : do many things for atoms
#  (5) valid_multiple_marginal : check the marginal's setting
#  (6) valid_multiple_weight   : weight of each measures
#  (7) valid_multiple_distance : distance between support and measure
#  (8) check_images            : for image barycenter



# (1) valid_weight --------------------------------------------------------
#' @keywords internal
#' @noRd
valid_weight <- function(wx, m, input.name, function.name){
  m = round(m)
  if ((length(wx)==0)&&(is.null(wx))){
    return(rep(1/m, m))
  } else {
    if (!is.vector(wx)){
      stop(paste0("* ",function.name," : input ",input.name," must be a vector."))
    }
    if (length(wx)!=m){
      stop(paste0("* ",function.name," : input ",input.name," should be of length ",m,"."))
    }
    if (any(wx<0)){
      stop(paste0("* ",function.name," : input ",input.name," should contain no negative values."))
    }
    return(wx/base::sum(wx))
  }
}

# (2) valid_distance ------------------------------------------------------
#' @keywords internal
#' @noRd
valid_distance <- function(D, input.name, function.name){
  cond1 = is.matrix(D)
  cond2 = all(D >= 0)
  cond3 = (!any(is.na(D)))
  cond4 = (!any(is.infinite(D)))
  
  if (cond1&&cond2&&cond3&&cond4){
    return(D)
  } else {
    stop(paste0("* ",function.name," : input ",input.name," should be a matrix of no negative/Inf/NaN values."))
  }
}

# (3) valid_matrixed ------------------------------------------------------
#' @keywords internal
#' @noRd
valid_matrixed <- function(data, fname){
  if (is.vector(data)){
    return(matrix(data, ncol=1))
  } else {
    if (is.matrix(data)){
      return(data)
    } else {
      dname = paste0("'",deparse(substitute(data)),"'")
      stop(paste0("* ",fname," : input ",dname," should be either a vector or a matrix."))
    }
  }
}


# (4) valid_multiple_measures ---------------------------------------------
#     Requirements : 1) list, 2) each is matrix, 3) all same dimension
#' @keywords internal
#' @noRd
valid_multiple_measures <- function(data, p, fname){
  dname = paste0("'",deparse(substitute(data)),"'")
  if (!is.list(data)){
    stop(paste0("* ",fname," : input ",dname," should be a list."))
  }
  K      = length(data)
  output = list()
  for (k in 1:K){
    tgt = data[[k]]
    if (is.vector(tgt)){
      output[[k]] = matrix(tgt, ncol=1)
    } else {
      if (is.matrix(tgt)){
        output[[k]] = tgt
      } else {
        stop(paste0("* ",fname," : ",k,"-th element of ",dname," should be either a vector or a matrix."))
      }
    }
  }
  for (k in 1:K){
    if (ncol(output[[k]])!=p){
      stop(paste0("* ",fname," : ",k,"-th element of ",dname," does not have the same dimension as the support atoms."))
    }
  }
  return(output)
}

# (5) valid_multiple_marginal : check the marginal's setting --------------
#' @keywords internal
#' @noRd
valid_multiple_marginal <- function(marginals, natoms, fname){
  dname = paste0("'",deparse(substitute(data)),"'")
  if ((length(marginals)==0)&&is.null(marginals)){
    K      = length(natoms)
    output = list()
    for (k in 1:K){
      nobjk       = natoms[k]
      output[[k]] = rep(1/nobjk, nobjk)
    }
    return(output)
  } else {
    K = length(natoms)
    if ((!is.list(marginals))||(length(marginals)!=K)){
      stop(paste0("* ",fname," : ",dname," should be a list of length ",K,"."))
    }
    output = list()
    for (k in 1:K){
      tgt = marginals[[k]]
      if ((length(tgt)!=natoms[k])||(any(tgt<0))){
        stop(paste0("* ",fname," : ",k,"-th element of ",dname," should be a nonnegative vector of length ",natoms[k],"."))
      }
      output[[k]] = tgt/sum(tgt)
    }
    return(output)
  }
}

# (6) valid_multiple_weight   : weight of each measures -------------------
#' @keywords internal
#' @noRd
valid_multiple_weight <- function(weight, K, fname){
  dname = paste0("'",deparse(substitute(weight)),"'")
  if ((length(weight)==0)&&is.null(weight)){
    return(rep(1/K, K))
  } else {
    if ((!is.vector(weight))||(length(weight)!=K)||(any(weight < 0))){
      stop(paste0("* ",fname," : ",dname," should be a length-",K," vector of nonnegative weights."))
    }
    return(weight/base::sum(weight))
  }
}

# (7) valid_multiple_distance ---------------------------------------------
#' @keywords internal
#' @noRd
valid_multiple_distance <- function(distances, fname){
  dname = paste0("'",deparse(substitute(distances)),"'")
  if (!is.list(distances)){
    stop(paste0("* ",fname," : ",dname," should be a length-",K," list."))
  }
  K = length(distances)
  output = list()
  for (k in 1:K){
    tgt = distances[[k]]
    if (is.vector(tgt)){
      output[[k]] = matrix(tgt, ncol=1)
    } else {
      if (!is.matrix(tgt)){
        stop(paste0("* ",fname," : ",k,"-th element of ",dname," should be a matrix of distances."))
      } else {
        output[[k]] = tgt
      }
    }
  }
  if (length(unlist(lapply(output,nrow)))!=1){
    stop(paste0("* ",fname," : not all elements of ",dname," have the same number of rows."))
  }
  return(output)
}

# (8) check_images --------------------------------------------------------
#' @keywords internal
#' @noRd
check_images <- function(images, fname){
  dname = paste0("'",deparse(substitute(images)),"'")
  if (!is.list(images)){
    stop(paste0("* ",fname," : ",dname," should be a list.")) 
  }
  if (any(unlist(lapply(images, is.matrix))==FALSE)){
    stop(paste0("* ",fname," : all elements of ",dname," should be a matrix."))
  }
  nnrows = unique(unlist(lapply(images, nrow)))
  nncols = unique(unlist(lapply(images, ncol)))
  if (length(nnrows)!=1){
    stop(paste0("* ",fname," : all image matrices should have same number of rows."))
  }
  if (length(nncols)!=1){
    stop(paste0("* ",fname," : all image matrices should have same number of columns."))
  }
  return(TRUE)
}