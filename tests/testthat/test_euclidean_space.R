context("Unit tests for Euclidean Space")

test_that("Tests instantiation of dimension", {
  rand <- ceiling(10*runif(1))
  euclidean.space <- EuclideanSpace$new(dimension = rand)
  expect_equal(euclidean.space$dimension,rand)
})

test_that("Tests RandomUniform",{
  rand <- ceiling(10*runif(1))
  euclidean.space <- EuclideanSpace$new(dimension = rand)
  point <- array(euclidean.space$RandomUniform())
  expect_true(length(point) == rand)
  expect_equivalent(c((point > -1) & (point < 1)), rep(TRUE,rand))
})

test_that("Tests point and belongs together", {
  rand <- ceiling(10*runif(1))
  euclidean.space <- EuclideanSpace$new(dimension = rand)
  point <- array(euclidean.space$RandomUniform())
  belongs <- euclidean.space$Belongs(point)
  expect_true(belongs)
})

test_that("Tests InnerProductMatrix",{
  rand <- ceiling(10*runif(1))
  euclidean.metric <- EuclideanMetric$new(dimension = rand)
  mat <- euclidean.metric$InnerProductMatrix()
  expected <- ToNdarray(diag(rand), to.ndim = 3)
  expect_equivalent(mat, expected)
})

