context("Unit tests for Vectorization")

test_that("Tests ToNdArray to 0 axis", {
  x <- 1:12 ; dim(x) <- c(3, 4)
  axis.0 <- ToNdarray(x, 3, 0)
  expect_equivalent(dim(axis.0), c(1, 3, 4))
})

test_that("Tests ToNdArray to 1 axis", {
  x <- 1:12 ; dim(x) <- c(3, 4)
  axis.1 <- ToNdarray(x, 3, 1)
  expect_equivalent(dim(axis.1), c(3, 1, 4))
})

test_that("Tests ToNdArray to 2 axis", {
  x <- 1:12 ; dim(x) <- c(3, 4)
  axis.2 <- ToNdarray(x, 3, 2)
  expect_equivalent(dim(axis.2), c(3, 4, 1))
})


