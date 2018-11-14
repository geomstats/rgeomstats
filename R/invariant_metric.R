InvariantMetric <- setRefClass("InvariantMetric",
  fields = c("group", "inner.product.mat.at.identity", "left.or.right"),
  methods = list(
    initialize = function(group, inner.product.mat.at.identity=NULL,
                          left.or.right="left") {
      if (is.null(inner.product.mat.at.identity)) {
        inner.product.mat.at.identity <- diag(1, group$dimension)
      }
      mat.shape = dim(inner.product.mat.at.identity)
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
    }
  )
)
