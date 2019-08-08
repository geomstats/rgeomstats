#' Riemannian and pseudo-Riemannian metrics.

RiemannianMetric <- setRefClass("RiemannianMetric",
  fields = c("object"),
  methods = list(

    initialize = function(dimension){
      stopifnot(dimension %% 1 == 0)
      stopifnot(dimension > 0)
      .self$dimension <- dimension
    },

    InnerProductMatrix = function(base.point = NULL){
      warning("The computation of the inner product matrix is not implemented.")
    },

    InnerProduct = function(tangent.vec.a, tangent.vec.b, base.point = NULL){
      tangent.vec.a <- ToNdarray(tangent.vec.a, to.ndim = 2)
      tangent.vec.b <- ToNdarray(tangent.vec.b, to.ndim = 2)
      n.tangent.vec.a <- dim(tangent.vec.a)[1]
      n.tangent.vec.b <- dim(tangent.vec.b)[1]

      inner.prod.mat <- .self$inner.product.matrix(base.point)
      inner.prod.mat <- ToNdarray(inner.prod.mat, to.ndim = 3)
      n.mats <- dim(inner.prod.mat)[1]

      n.inner.prod <- max(n.tangent.vec.a, n.tangent.vec.b)
      n.inner.prod <- max(n.inner.prod, n.mats)

      n.tiles.a <- n.inner.prod / n.tangent.vec.a
      tangent.vec.a <- rep(tangent.vec.a, n.tiles.a)

      n.tiles.b <- n.inner.prod / n.tangent.vec.b
      tangent.vec.b <- rep(tangent.vec.b, n.tiles.b)

      n.tiles.mat <- n.inner.prod / n.mats
      inner.prod.mat <- rep(inner.prod.mat, n.tiles.mat)
      inner.prod.mat <- ToNdarray(inner.prod.mat, axis = 3)

      aux <- tangent.vec.a %*% inner.prod.mat
      inner.prod <- aux * tangent.vec.b
      inner.prod <- ToNdarray(inner.prod, to.ndim = 2, axis = 1)

      stopifnot(length(dim(inner.prod) == 2))
      return(inner.prod)
    },

    SquaredNorm = function(vector, base.point= NULL){
      sq.norm <- .self$InnerProduct(vector, vector, base.point)
      return(sq.norm)
    },

    Norm = function(vector, base.point=NULL){
      sq.norm <- .self$SquaredNorm(vector, base.point)
      norm <- sqrt(sq.norm)
      return(norm)
    }
  )
)
