# Unit tests for Euclidean Space

setUp <- function(){
  rand <<- ceiling(10*runif(1))
  euclidean.space <<- EuclideanSpace$new(dimension = rand)
}

test_that("instantiating random dimension 1-10", {
  rand <<- ceiling(10*runif(1))
  euclidean.space <<- EuclideanSpace$new(dimension = rand)

  expect_equal(euclidean.space$dimension,rand)
})


test_that("RandomUniform and Belongs together", {
  point <- array(euclidean.space$RandomUniform())
  belongs <- euclidean.space$Belongs(point)
  expect_true(belongs)
})
