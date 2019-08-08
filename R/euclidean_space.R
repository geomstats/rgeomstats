#' Euclidean space.

#' @include manifold.R
EuclideanSpace <- setRefClass("EuclideanSpace",
                              contains = "Manifold",
                              fields = "dimension",
                              methods = list(
                                initialize = function(dimension){
                                  stopifnot(dimension %% 1 == 0)
                                  stopifnot(dimension > 0)
                                  .self$dimension <- dimension
                                },

                                Belongs = function(point){
                                  "Evaluate if a point belongs to the Euclidean space."
                                  point <- ToNdarray(point, to.ndim = 2, axis = 0)
                                  n.points <- dim(point)[1]
                                  points.dim <- dim(point)[2]
                                  belongs <- points.dim == .self$dimension
                                  belongs <- ToNdarray(belongs, to.ndim = 1, axis = 0)
                                  belongs <- ToNdarray(belongs, to.ndim = 2, axis = 1)
                                  belongs <- rep(belongs, n.points)
                                  return(belongs)
                                },

                                RandomUniform = function(n_samples = 1){
                                  "Sample in the Euclidean space with the uniform distribution."
                                  point <- replicate(dimension, (runif(n_samples) - 0.5) * 2)
                                  return(point)
                                }

                              )
)


EuclideanMetric <- setRefClass("EuclideanMetric",
                               fields = "dimension",
                               methods = list(
                                 initialize = function(dimension){
                                   stopifnot(dimension %% 1 == 0)
                                   stopifnot(dimension > 0)
                                   .self$dimension <- dimension
                                 },

                                 InnerProductMatrix = function(base.point=None){
                                   "Inner product matrix, independent of the base point."
                                   mat <- diag(.self$dimension)
                                   mat <- ToNdarray(mat, to.ndim = 3)
                                   return(mat)
                                 },


                                 Exp = function(tangent.vec, base.point){
                                   "The Riemannian exponential is the addition in the Euclidean space."
                                   tangent.vec <- ToNdarray(tangent.vec, to.ndim = 2)
                                   base.point <- ToNdarray(base.point, to.ndim = 2)
                                   riemexp <- base.point + tangent.vec
                                   return(riemexp)
                                 },

                                 Log = function(point, base.point){
                                   "The Riemannian logarithm is the subtraction in the Euclidean space."
                                   point <- ToNdarray(base.point, to.ndim = 2)
                                   base.point <- ToNdarray(base.point, to.ndim = 2)
                                   riemlog <- point - base.point
                                   return(riemlog)
                                 },

                                 Mean = function(points, weights=None){
                                   "The Frechet mean of (weighted) points computed with the
                                   Euclidean metric is the weighted average of the points
                                   in the Euclidean space."
                                   riemmean <- weighted.mean(points, weights = weights)
                                   riemmean <- ToNdarray(mean, to.ndim = 2)
                                   return(riemmean)
                                 }
                               )

)
