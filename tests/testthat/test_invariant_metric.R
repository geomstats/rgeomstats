context("Unit tests for Invariant Metric")

test_that("Tests instantiation of Invariant Metric", {
  so3 <- SpecialOrthogonalGroup$new(n = 3)
  invariant.metric <- InvariantMetric$new(group = so3)
  expect_equivalent(invariant.metric$inner.product.mat.at.identity, diag(1, 3))
  expect_equivalent(invariant.metric$left.or.right, "left")
})

