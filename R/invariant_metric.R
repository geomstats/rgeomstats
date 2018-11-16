InvariantMetric <- setRefClass("InvariantMetric",
  fields = c("group", "inner.product.mat.at.identity", "left.or.right"),
  contains = "RiemannianMetric",
  methods = list(
    initialize = function(group, inner.product.mat.at.identity=NULL,
                          left.or.right = "left") {
      if (is.null(inner.product.mat.at.identity)) {
        inner.product.mat.at.identity <- diag(1, group$dimension)
      }
      mat.shape <- dim(inner.product.mat.at.identity)
      stopifnot(left.or.right == "left" || "right")
      eigenvalues <- eigen(inner.product.mat.at.identity)
      eigenvalues <- eigenvalues$values
      n.pos.eigval <- sum(eigenvalues > 0)
      n.neg.eigval <- sum(eigenvalues < 0)
      n.null.eigval <- sum(eigenvalues < .Machine$double.eps ^ 0.5)
      .self$group <- group

      if (is.null(inner.product.mat.at.identity)) {
        inner.product.mat.at.identity <- diag(1, .self$group$dimension)
      }
      .self$inner.product.mat.at.identity <- inner.product.mat.at.identity
      .self$left.or.right <- left.or.right
    },

    InnerProductAtIdentity = function(tangent.vec.a, tangent.vec.b){
      tangent.vec.a <- ToNdarray(tangent.vec.a, to.ndim = 2)
      tangent.vec.b <- ToNdarray(tangent.vec.b, to.ndim = 2)

      inner.prod <- (tangent.vec.a %*% .self$inner.product.mat.at.identity) %*% t(tangent.vec.b)
      inner.prod <- ToNdarray(inner.prod, to.ndim = 2, axis = 1)

      return(inner.prod)
    },

    InnerProduct = function(tangent.vec.a, tangent.vec.b, base.point = NULL){
    if (is.null(base.point)) {
      return(.self$InnerProductAtIdentity(tangent.vec.a, tangent.vec.b))
    } else{
      return(callSuper(InnerProduct(tangent.vec.a, tangent.vec.b, base.point)))
    }

    if (.self$left.or.right == "right") {
      warning('inner.product not implemented for right invariant metrics.')
    }
    jacobian <- .self$group.JacobianTranslation(base.point)
    inv.jacobian <- solve(jacobian)
    tangent.vec.a.at.id <- inv.jacobian %*% tangent.vec.a
    tangent.vec.b.at.id <- inv.jacobian %*% tangent.vec.b
    inner.prod <- .self$InnerProductAtIdentity(tangent.vec.a.at.id,
                                                tangent.vec.b.at.id)
    return(inner.prod)
    }
  )
)
