test_that("Unit tests for Euclidean Space", {
  rand <<- ceiling(10*runif(1))
  euclidean.space <<- EuclideanSpace$new(dimension = rand)
  expect_equal(euclidean.space$dimension,rand)
})


test_that("Unit tests for Euclidean Space", {
  point <- array(euclidean.space$RandomUniform())
  belongs <- euclidean.space$Belongs(point)
  expect_true(belongs)
})
