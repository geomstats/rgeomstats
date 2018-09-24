# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#roxygen2 documentation will be run with roxygen2::roxygenise()
#roxygen2 generates mans with documentation embedded next to code
#documentation found in man folder and functions by ?<function>

#' Function name
#'
#' more detailed description
#'
#' @param x numeric vector
#'
#' @return hello, world
#'
#'
#' @examples
#' hello()

hello <- function() {
  print("Hello, world!")
}
