# Manifold, i.e. a topological space that locally resembles
# Euclidean space near each point.

library(reticulate)
gs <- import_from_path("geomstats", path = ".")

Manifold <- setRefClass("Manifold",
                        fields = "dimension",
                        methods = list(
                          initialize = function(dimension){
                            stopifnot(dimension %% 1 == 0)
                            stopifnot(dimension > 0)
                            .self$dimension <- dimension
                          }
                        )
)
