
#Install reticulate first time
install.packages("reticulate")


library(reticulate)
py_install("pandas")

# install numpy
conda_install("r-reticulate", "numpy")

# import numpy (it will be automatically discovered in "r-reticulate")
np <- import("numpy")


#' flip
#'
#' @param m Input array
#' @param axis None or int or tuple of ints, optional
#'
#' @return array_like A view of m with the entries of axis reversed. Since a view is returned, this operation is done in constant time.
flip = function(...){
  return(np$flip(...))
}


#' amax
#'
#' @param axis None or int or tuple of ints, optional
#' @param out ndarray, optional
#' @param keepdims bool, optional
#' @param initial scalar, optional
#'
#' @return amax Maximum of a. If axis is None, the result is a scalar value.
amax = function(...){
  return(np$amax(...))
}




#' arctan2
#'
#' @param x1 array_like, real-valued
#' @param x2 array_like, real-valued
#' @param out ndarray, None, or tuple of ndarray and None, optional
#' @param where array_like, optional
#'
#' @return angle Array of angles in radians, in the range [-pi, pi]. This is a scalar if both x1 and x2 are scalars.
arctan2 = function(...){
  return(np$arctan2(...))
}

