# Unit tests for Euclidean Space

dimension=3;
length(euclidean.space$RandomUniform()) ==

SetUp <- function(){
  rand <<- ceiling(10*runif(1))
  euclidean.space <<- EuclideanSpace$new(dimension = rand)
}

TestDimension <- function(){
  result <<- euclidean.space$dimension
  expected <<- rand
  stopifnot(result == expected)
}

TestBelongs <- function(){
  euclidean.space$Belongs(c(runif(euclidean.space$dimension)))
}
