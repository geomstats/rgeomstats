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

                                          if (.self$n == 3){
                                            # applys norm to each row
                                            angle <- apply(regularized.point, 1, function(x){sqrt(sum(x^2))})
                                            mask.0 <- (angle < .Machine$double.eps ^ 0.5)
                                            mask.not.0 <- !mask.0
                                            mask.pi <- ((angle - pi) < .Machine$double.eps ^ 0.5)

                                            k <- floor(angle / (2 * pi) + .5)

                                            norms.ratio <- array(0,c(n.points,1))

                                            # This avoids division by 0.
                                            angle <- angle + mask.0 * 1.

                                            norms.ratio <- norms.ratio + mask.not.0 * (1 - 2 * pi * k / angle)
                                            norms.ratio <- norms.ratio + mask.0 * 1.
                                            norms.ratio <- norms.ratio + mask.pi * (pi / angle
                                                                                          - (1 - 2 * pi * k / angle))

                                            regularized.point <- apply(regularized.point, 2, function(x){x * norms.ratio})
                                            if (is.null(dim(regularized.point))) {
                                              regularized.point <- ToNdarray(array(regularized.point), to.ndim = 2)
                                            }
                                          }
                                          stopifnot(length(dim(regularized.point)) == 2)
                                          return(regularized.point)
                                        }
                                      )
)
