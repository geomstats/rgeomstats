# euclidean_space

library(reticulate)
gs <- import_from_path("geomstats", path = ".")

EuclideanSpace <- setRefClass("EuclideanSpace",
                              contains = "Manifold",
                              methods = list(
                                
                                initialize <- function(self, dimension){
                                  stopifnot(dimension %% 1 == 0, dimension > 0)
                                  self$dimension <- dimension
                                },
                                
                                Belongs <- function(self, point){
                                  "Evaluate if a point belongs to the Euclidean space."
                                  point <- gs$backend$to_ndarray(point, to_ndim = 2)
                                  n_points = point$shape
                                  points_dim = point$shape
                                  belongs = point_dim == self$dimension
                                  belongs = gs$backend$to_ndarray(belongs, to_ndim = 1)
                                  belongs = gs$backend$to_ndarray(belongs, to_ndim = 2, axis = 1)
                                  belongs = gs$backend$tile(belongs, c(n_points, 1))
                                  return(belongs)
                                },
                                
                                RandomUniform <- function(self, n_samples = 1){
                                  "Sample in the Euclidean space with the uniform distribution."
                                  size <- c(n_samples, self$dimension)
                                  point <- (gs$backend$random$rand(size) - 0.5) * 2
                                  return(point)
                                }
                                
                              )
)


EuclideanMetric <- setRefClass("EuclideanMetric",
                               methods = list(
                                 
                                 Initialize <- function(self, dimension){
                                   stopifnot(dimension %% 1 == 0, dimension > 0)
                                 },
                                 
                                 InnerProductMatrix <- function(self, base_point=None){
                                   "Inner product matrix, independent of the base point."
                                   
                                   mat <- diag(self$dimension)
                                   mat <- gs$backend$to_ndarray(mat, to_ndim = 3)
                                   return(mat)
                                 },
                                 
                                 
                                 Exp <- function(self, tangent_vec, base_point){
                                   "The Riemannian exponential is the addition in the Euclidean space."
                                   tangent_vec <- gs$backend$to_ndarray(tangent_vec, to_ndim = 2)
                                   base_point <- gs$backend$to_ndarray(base_point, to_ndim = 2)
                                   Riemexp <- base_point + tangent_vec
                                   return(Riemexp)
                                 },
                                 
                                 Log <- function(self, point, base_point){
                                   "The Riemannian logarithm is the subtraction in the Euclidean space."
                                   point <- gs$backend$to_ndarray(base_point, to_ndim = 2)
                                   base_point <- gs$backend$to_ndarray(base_point, to_ndim = 2)
                                   Riemlog <- point - base_point
                                   return(Riemlog)
                                 },
                                 
                                 Mean <- function(self, points, weights=None){
                                   "The Frechet mean of (weighted) points computed with the
                                   Euclidean metric is the weighted average of the points
                                   in the Euclidean space."
                                   Riemmean <- gs$backend$average(points, axis = 0, weights = weights)
                                   Riemmean <- gs$backend$to_ndarray(mean, to_ndim = 2)
                                   return(Riemmean)
                                 }
                               )
                               
)