#Z-norm is frequently used in other function so it is defined as a separate function outside
#for frequent re-use within the package

#' z-normalise data
#'
#' Function applicable on any numeric distribution of data and z-normalize them.
#'
#' @param data Numeric distribution; may be stored in an object of type vector, matrix or data frame
#'
#' @return z-normalized distribution of input
#'
#' @examples
#' random_distribution <- runif(n = 50, min = 1, max = 10)
#' znorm_distribution <- znorm(random_distribution)
#'
#' @export
znorm <- function (data) {
  data.mean <- mean(data)
  data.dev <- sd(data)
  (data - data.mean)/data.dev
}
