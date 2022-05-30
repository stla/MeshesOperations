isAtomicVector <- function(x){
  is.atomic(x) && is.vector(x)
}

isPositiveNumber <- function(x){
  is.numeric(x) && length(x) == 1L && x > 0 && !is.na(x)
}

isNonNegativeNumber <- function(x){
  is.numeric(x) && length(x) == 1L && x >= 0 && !is.na(x)
}

isPositiveInteger <- function(x){
  is.numeric(x) && length(x) == 1L && !is.na(x) && floor(x) == x
}

isStrictPositiveInteger <- function(x){
  isPositiveInteger(x) && x > 0
}

isBoolean <- function(x){
  is.logical(x) && length(x) == 1L && !is.na(x)
}

#' @importFrom gmp `%*%`
#' @noRd
`%^%` <- function(A, n){
  Reduce(gmp::`%*%`, replicate(n, A, simplify = FALSE))
}
