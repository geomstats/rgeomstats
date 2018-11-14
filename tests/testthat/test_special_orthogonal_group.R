context("Unit tests for Special Orthogonal Group")

elements <- list(
  with.angle.0 = array(c(0,0,0)),
  with.angle.close.0 = 1e-10 * array(c(1, -1, 1)),
  with.angle.close.pi.low = ((pi - 1e-9) / sqrt(2) * array(c(0, 1, -1))),
  with.angle.pi = pi / sqrt(3) * array(c(1, 1, -1)),
  with.angle.close.pi.high = ((pi + 1e-9) / sqrt(3) * array(c(-1, 1, -1))),
  with.angle.in.pi.2pi = ((pi + 0.3) / sqrt(5) * array(c(-2, 1, 0))),
  with.angle.close.2pi.low = ((2 * pi - 1e-9) / sqrt(6) * array(c(2, 1, -1))),
  with.angle.2pi = 2 * pi / sqrt(3) * array(c(1, 1, -1)),
  with.angle.close.2pi.high = ((2 * pi + 1e-9) / sqrt(2) * array(c(1, 0, -1)))
)

test_that("Tests instantiation of dimension", {
  rand <- ceiling(10 * runif(1))
  expected.dimension <- (rand * (rand - 1)) / 2
  so3 <- SpecialOrthogonalGroup$new(n = rand)
  expect_equal(so3$n, rand)
  expect_equal(so3$dimension, expected.dimension)
})

test_that("Tests Regularize", {
  so3 <- SpecialOrthogonalGroup$new(n = 3)
  expect_equivalent(so3$Regularize(point = elements$with.angle.close.0),
                    ToNdarray(elements$with.angle.close.0, to.ndim = 2))
  expect_equivalent(so3$Regularize(point = elements$with.angle.close.pi.low),
                    ToNdarray(elements$with.angle.close.pi.low, to.ndim = 2))
  expect_equivalent(so3$Regularize(point = elements$with.angle.pi),
                    ToNdarray(elements$with.angle.pi, to.ndim = 2))
  expect_equivalent(so3$Regularize(point = elements$with.angle.close.pi.high),
                    ToNdarray(elements$with.angle.close.pi.high, to.ndim = 2))
  expect_equivalent(so3$Regularize(point = elements$with.angle.2pi),
                    ToNdarray(array(c(0,0,0)), to.ndim = 2))


  # for angle between pi and 2pi
  angle <- sqrt(sum(elements$with.angle.in.pi.2pi ^ 2))
  new.angle <- pi - (angle - pi)
  expected <- -(new.angle / angle) * elements$with.angle.in.pi.2pi

  expect_equivalent(so3$Regularize(point = elements$with.angle.in.pi.2pi),
                    ToNdarray(expected,to.ndim = 2))

  # for angle 2pi low
  angle <- sqrt(sum(elements$with.angle.close.2pi.low ^ 2))
  new.angle <- pi - (angle - pi)
  expected <- -(new.angle / angle) * elements$with.angle.close.2pi.low

  expect_equivalent(so3$Regularize(point = elements$with.angle.close.2pi.low),
                    ToNdarray(expected,to.ndim = 2))

  # for angle 2pi high
  angle <- sqrt(sum(elements$with.angle.close.2pi.high^2))
  new.angle <- angle - 2 * pi
  expected <- new.angle * elements$with.angle.close.2pi.high / angle
  expect_equivalent(so3$Regularize(point = elements$with.angle.close.2pi.high),
                    ToNdarray(expected,to.ndim = 2))


})

test_that("Tests Skew Matrix From Vector", {
  so3 <- SpecialOrthogonalGroup$new(n = 3)
  vec <- ToNdarray(array(runif(3)), to.ndim = 2)
  result <- array(so3$SkewMatrixFromVector(vec), dim = c(3, 3))
  expect_equivalent(vec %*% result, ToNdarray(array(c(0, 0, 0)), to.ndim = 2))
})


test_that("Tests Jacobian Translation", {
  so3 <- SpecialOrthogonalGroup$new(n = 3)
  for (array.type in elements) {
    point <- array.type
    jacobian <- so3$JacobianTranslation(ToNdarray(array(point), to.ndim = 2))
    result <- det(array(jacobian, dim = c(3, 3)))
    point <- so3$Regularize(point)
    angle <- sqrt(sum(point ^ 2))
    if (identical(array.type, elements$with.angle.0) ||
        identical(array.type, elements$with.angle.close.0) ||
        identical(array.type, elements$with.angle.2pi) ||
        identical(array.type, elements$with.angle.2pi.high)) {
      expected <- 1 + angle ^ 2 / 12 + angle ^ 4 / 240
    } else {
      expected <- angle ^ 2 / (4 * sin(angle / 2) ^ 2)
    }
    expect_equal(result, expected)
  }
})

test_that("Tests Vector From Skew Matrix", {
  so3 <- SpecialOrthogonalGroup$new(n = 3)
  vec <- ToNdarray(array(runif(3)), to.ndim = 2)
  skew.mat <- so3$SkewMatrixFromVector(vec)
  vec.from.skew.matrix <- so3$VectorFromSkewMatrix(skew.mat = skew.mat)
  expect_equivalent(vec, vec.from.skew.matrix)
})
