#' Utils to vectorize rgeomstats code
ToNdarray <- function(x, to.ndim, axis = 0) {
  CurrentNdim <- length(dim(x))
  if (CurrentNdim == to.ndim - 1) {
    if (axis == 0) {
      NewDim <- c(1, dim(x))
      NewX <- array(x, dim = NewDim)
    } else if (axis == CurrentNdim) {
      NewDim <- c(dim(x), 1)
      NewX <- array(x, dim = NewDim)
    } else {
      NewDim <- c(dim(x)[1:axis], 1, dim(x)[(axis + 1):length(dim(x))])
      NewX <- array(x, dim = NewDim)
    }
  } else if (CurrentNdim == to.ndim) {
    NewX <- x
  } else {
    warning(cat("We can increase ndim only by 1. Here: ndim of x is:", length(dim(x))))
  }
  return(NewX)
}
