context("Unit tests for Invariant Metric")

test_that("Tests instantiation of Invariant Metric", {
  so3 <- SpecialOrthogonalGroup$new(n = 3)
  invariant.metric <- InvariantMetric$new(group = so3)
  expect_equivalent(invariant.metric$inner.product.mat.at.identity, diag(1, 3))
  expect_equivalent(invariant.metric$left.or.right, "left")
})

test_that("Tests inner product and norm",{
  so3 <- SpecialOrthogonalGroup$new(n = 3)
  im <- InvariantMetric$new(so3)
  tangent.vec.a <- array(c(1, 2, 3))
  tangent.vec.b <- array(c(3, 2, 1))

  result <- im$InnerProduct(tangent.vec.a, tangent.vec.b)
  expected <- 10
  expect_equal(result[1], expected)

  result <- im$Norm(tangent.vec.a)
  expected <- sqrt(14)
  expect_equal(result[1], expected)
})
