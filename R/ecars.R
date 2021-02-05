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
ecars <- function(X, Y, alpha=NULL){
  if(is.null(alpha))
  {
    alpha <- median(alpha_estimate(X, Y))
  }
  corpcor::crossprod.powcor.shrink(X, Y, alpha)
}
