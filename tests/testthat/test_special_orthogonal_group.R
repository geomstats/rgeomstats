context("Unit tests for Special Orthogonal Group")

test_that("Tests instantiation of dimension", {
  rand <- ceiling(10*runif(1))
  expected.dimension <- (rand * (rand - 1)) / 2
  special.orthogonal.group <- SpecialOrthogonalGroup$new(n = rand)
  expect_equal(special.orthogonal.group$n, rand)
  expect_equal(special.orthogonal.group$dimension, expected.dimension)
})

