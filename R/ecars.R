#' Calculate importance
#'
#' Assign importance value to each variable in the
#' regression model.
#'
#' @param X A nxp data matrix
#' @param Y A vector
#' @param alpha estimated alpha
#'
#' @return A vector of each variable's importance value
#' @export
#'
#' @examples
#' ecars(X)
ecars <- function(X, Y, alpha){
  corpcor::crossprod.powcor.shrink(X, Y, alpha)
}
