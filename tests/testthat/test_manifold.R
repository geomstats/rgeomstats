context("Unit tests for Manifold")

test_that("Tests instantiation of dimension",{
  rand <- ceiling(10*runif(1))
  manifold <- Manifold$new(dimension = rand)
  result <- manifold$dimension
  expected <- rand
  expect_equal(result, expected)
})
