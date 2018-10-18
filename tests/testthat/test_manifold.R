
test_that("Unit tests for Manifold",{
  rand <- ceiling(10*runif(1))
  manifold <- Manifold$new(dimension = rand)
  result <- manifold$dimension
  expected <- rand
  expect_equal(result, expected)
})
