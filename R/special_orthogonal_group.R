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
                                        }
                                      )
)
