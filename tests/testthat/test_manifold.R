# Unit tests for Manifold.

library(reticulate)
gs <- import_from_path("geomstats", path = ".")

  rand <<- ceiling(10*runif(1))
  manifold <<- Manifold$new(dimension = rand)



test_that("testing dimension",{
  TestDimension <- function(){
    result <<- manifold$dimension
    expected <<- rand
    test.dimension = result == expected
    return(test.dimension)
  }
  expect_true(TestDimension())
})

expect_true()

TestBelongs <- function(){
  point <- euclidean.space$RandomUniform()
  belongs <- euclidean.space$Belongs(point)
}
