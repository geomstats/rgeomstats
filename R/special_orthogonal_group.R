# The special orthogonal group SO(n),
# i.e. the Lie group of rotations in n dimensions.


kATOL <- 1e-5

TAYLOR.COEFFS.1.AT.0 <- c(1., 0.,
                          - 1. / 12., 0.,
                          - 1. / 720., 0.,
                          - 1. / 30240., 0.)
TAYLOR.COEFFS.2.AT.0 <- c(1. / 12., 0.,
                          1. / 720., 0.,
                          1. / 30240., 0.,
                          1. / 1209600., 0.)
TAYLOR.COEFFS.1.AT.PI <- c(0., - pi / 4.,
                           - 1. / 4., - pi / 48.,
                           - 1. / 48., - pi / 480.,
                           - 1. / 480.)

SpecialOrthogonalGroup <- setRefClass("SpecialOrthogonalGroup",
  fields = c("n", "dimension"),
  methods = list(
    initialize = function(n){
      stopifnot(n %% 1 == 0)
      stopifnot(n > 0)
      .self$n <- n
      .self$dimension <- (n * (n - 1)) / 2
    },

    Belongs = function(point){
      "Evaluate if a point belongs to SO(n)."
      # for point type vector

      point <- ToNdarray(point, to.ndim = 2)
      vec.dim <-  dim(point)[2]

      return(vec.dim == .self$dimension)
    },

    Regularize = function(point){
      "In 3D, regularize the norm of the rotation vector,
      to be between 0 and pi, following the axis-angle
      representation's convention.
      If the angle angle is between pi and 2pi,
      the function computes its complementary in 2pi and
      inverts the direction of the rotation axis."
      point <- ToNdarray(point, to.ndim = 2)
      n.points <- dim(point)[1]
      regularized.point <- point

      if (.self$n == 3) {
        angle <- apply(regularized.point, 1, function(x){sqrt(sum(x ^ 2))})
        mask.0 <- (angle < .Machine$double.eps ^ 0.5)
        mask.not.0 <- !mask.0
        mask.pi <- ((angle - pi) < .Machine$double.eps ^ 0.5)

        k <- floor(angle / (2 * pi) + .5)

        norms.ratio <- array(0, c(n.points, 1))

        # This avoids division by 0.
        angle <- angle + mask.0 * 1.

        norms.ratio[mask.not.0] <- 1 - 2 * pi * k / angle
        norms.ratio[mask.0] <- 1
        norms.ratio[mask.pi] <-  (pi / angle[mask.pi])

        regularized.point <- apply(regularized.point, 2, function(x){x * norms.ratio})
        if (is.null(dim(regularized.point))) {
          regularized.point <- ToNdarray(array(regularized.point), to.ndim = 2)
        }
      }
      stopifnot(length(dim(regularized.point)) == 2)
      return(regularized.point)
    },

    SkewMatrixFromVector = function(vec){
      "In 3D, compute the skew-symmetric matrix,
      known as the cross-product of a vector,
      associated to the vector vec."
      vec <- ToNdarray(vec, to.ndim = 2)
      n.vecs <- dim(vec)[1]
      vec.dim <- dim(vec)[2]
      if (.self$n == 3) {
        levi.civita.symbol <- array(data = c(
          c(0, 0, 0),
          c(0, 0, 1),
          c(0, -1, 0),
          c(0, 0, -1),
          c(0, 0, 0),
          c(1, 0, 0),
          c(0, 1, 0),
          c(-1, 0, 0),
          c(0, 0, 0)
        ), dim = c(3, 3, 3)
        )

        basis.vec.1 <- ToNdarray(array(c(1, 0, 0) * n.vecs), to.ndim = 2)
        basis.vec.2 <- ToNdarray(array(c(0, 1, 0) * n.vecs), to.ndim = 2)
        basis.vec.3 <- ToNdarray(array(c(0, 0, 1) * n.vecs), to.ndim = 2)

        ApplyLeviCivitaSymbol <- function(vec) {
          array(cbind(
          vec %*% levi.civita.symbol[1, , ],
          vec %*% levi.civita.symbol[2, , ],
          vec %*% levi.civita.symbol[3, , ]),
          dim = c(3, 3))
        }

        cross.prod.1 <- vec %*% ApplyLeviCivitaSymbol(basis.vec.1)
        cross.prod.2 <- vec %*% ApplyLeviCivitaSymbol(basis.vec.2)
        cross.prod.3 <- vec %*% ApplyLeviCivitaSymbol(basis.vec.3)

        skew.mat <- array(c(cross.prod.1, cross.prod.2, cross.prod.3), dim = c(1, 3, 3))
      }
      stopifnot(length(dim(skew.mat)) == 3)
      return(skew.mat)
    },

    JacobianTranslation = function(point, left.or.right="left"){
      "Compute the jacobian matrix of the differential
      of the left/right translations from the identity to point in SO(n)."

      stopifnot(left.or.right == "left" || "right")

      if (.self$n == 3) {
        point <- .self$Regularize(point)
        n.points <- dim(point)[1]
        angle <- sqrt(sum(point ^ 2))
        angle <- ToNdarray((array(angle)), to.ndim = 2, axis = 1)

        coef.1 <- array(0, dim = c(n.points, 1))
        coef.2 <- array(0, dim = c(n.points, 1))

        mask.0 <- (angle < .Machine$double.eps ^ 0.5)

        coef.1[mask.0] <- (
          TAYLOR.COEFFS.1.AT.0[1]
          + TAYLOR.COEFFS.1.AT.0[3] * angle[mask.0] ^ 2
          + TAYLOR.COEFFS.1.AT.0[5] * angle[mask.0] ^ 4
          + TAYLOR.COEFFS.1.AT.0[7] * angle[mask.0] ^ 6)

        coef.2[mask.0] <- (
          TAYLOR.COEFFS.2.AT.0[1]
          + TAYLOR.COEFFS.2.AT.0[3] * angle ^ 2
          + TAYLOR.COEFFS.2.AT.0[5] * angle ^ 4
          + TAYLOR.COEFFS.2.AT.0[7] * angle ^ 6)

        mask.pi <- (abs(angle - pi) < .Machine$double.eps ^ 0.5)

        delta.angle <- angle - pi

        coef.1[mask.pi] <- (
          TAYLOR.COEFFS.1.AT.PI[2] * delta.angle
          + TAYLOR.COEFFS.1.AT.PI[3] * delta.angle ^ 2
          + TAYLOR.COEFFS.1.AT.PI[4] * delta.angle ^ 3
          + TAYLOR.COEFFS.1.AT.PI[5] * delta.angle ^ 4
          + TAYLOR.COEFFS.1.AT.PI[6] * delta.angle ^ 5
          + TAYLOR.COEFFS.1.AT.PI[7] * delta.angle ^ 6)
        coef.2[mask.pi] <- (1 - coef.1[mask.pi]) / angle[mask.pi] ^ 2

        mask.else <- !mask.0 & !mask.pi
        angle <- angle + mask.pi
        coef.1[mask.else] <- ((angle[mask.else] / 2)
                              / tan(angle[mask.else] / 2))
        coef.2[mask.else] <- ((1 - coef.1[mask.else])
                              / angle[mask.else] ^ 2)
        jacobian <- array(0, dim = c(n.points, .self$dimension, .self$dimension))

        for (i in range(n.points)[1]:range(n.points)[2]) {
          sign <- -1
          if (left.or.right == "left") {
            sign <- 1
          }

          jacobian[i, , ] <- (
            coef.1[i] * diag(1, .self$dimension, .self$dimension)
            + coef.2[i] * outer(point[i, ], point[i, ])
            + sign * .self$SkewMatrixFromVector(ToNdarray(array(point[1, ]), to.ndim = 2))[i, , ] / 2)
        }
      }
      stopifnot(length(dim(jacobian)) == 3)
      return(jacobian)
    },

    VectorFromSkewMatrix = function(skew.mat){
      "In 3D, compute the vector defining the cross product
      associated to the skew-symmetric matrix skew mat.

      In nD, fill a vector by reading the values
      of the upper triangle of skew_mat."

      skew.mat <- ToNdarray(skew.mat, to.ndim = 3)
      n.skew.mats <- dim(skew.mat)[1]
      mat.dim.1 <- dim(skew.mat)[2]
      mat.dim.2 <- dim(skew.mat)[3]

      stopifnot(mat.dim.1 == mat.dim.2)
      stopifnot(mat.dim.2 == .self$n)

      vec.dim <- .self$dimension
      vec <- array(0, dim = c(n.skew.mats, vec.dim))

      if (.self$n == 3) {
        vec.1 <- ToNdarray(array(skew.mat[ , 3, 2]), to.ndim = 2, axis = 1)
        vec.2 <- ToNdarray(array(skew.mat[ , 1, 3]), to.ndim = 2, axis = 1)
        vec.3 <- ToNdarray(array(skew.mat[ , 2, 1]), to.ndim = 2, axis = 1)
        vec <- rbind(c(vec.1, vec.2, vec.3))
      }
      stopifnot(length(dim(vec)) == 2)
      return(vec)
    }
  )
)
