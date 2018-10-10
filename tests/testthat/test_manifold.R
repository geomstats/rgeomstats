# Unit tests for Manifold.

library(reticulate)
gs <- import_from_path("geomstats", path = ".")


manifold <- Manifold$new(dimension = 3)
manifold$dimension == 3

SetUp <- function(){
  rand <<- ceiling(10*runif(1))
  manifold <<- Manifold$new(dimension = rand)
}

TestDimension <- function(){
  result <<- manifold$dimension
  expected <<- rand
  stopifnot(result == expected)
}
