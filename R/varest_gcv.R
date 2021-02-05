#' Title Estimate error variance
#'
#' Use LASSO to estimate the error variance in the data
#' @param x An nxp data matrix
#' @param y A vector of response variable
#'
#' @return the estimated error variance
#' @export
#'
#' @examples
#' Sigma <- matrix(0.1, 500, 500)
#' X <- MASS::mvrnorm(500, rep(0, 500), Sigma)
#' ind <- 1:50
#' Y <- apply(X[, ind],2,sum) + rnorm(500,0, sqrt(5))
#' varest_gcv(X, Y)
varest_gcv <- function (x, y) {
  n = nrow(x)
  p = ncol(x)
  half = ceiling(n/2)
  x1 = x[1:half, ]
  y1 = y[1:half]
  x2 = x[-c(1:half), ]
  y2 = y[-c(1:half)]
  fit = glmnet::glmnet(x1, y1, family = "gaussian")
  ind = gcv_lambda_min(x1,y1,return="index")
  if (ind == 1){
    sigmahat1 = var(y2)
  }else {
    ind1 = which(abs(fit$beta[, ind]) > 0)
    if (length(ind1) > half/2) {
      betasort = sort(abs(fit$beta[ind1, ind]), decreasing = TRUE,
                      index.return = T)
      ind1 = ind1[betasort$ix[1:ceiling(half/2)]]
    }
    xm12 = x2[, ind1]
    sigmahat1 =  sum((y2 - lm(y2~xm12)$fitted.values)^2)/(half -
                                                            length(ind1)-1)
  }
  fit = glmnet::glmnet(x2, y2, family = "gaussian")
  ind = gcv_lambda_min(x2,y2,return="index")

  if (ind == 1) {
    sigmahat2 = var(y2)
  }else {
    ind1 = which(abs(fit$beta[, ind]) > 0)
    if (length(ind1) > half/2) {
      betasort = sort(abs(fit$beta[ind1, ind]), decreasing = TRUE,
                      index.return = T)
      ind1 = ind1[betasort$ix[1:ceiling(half/2)]]
    }
    xm21 = x1[, ind1]
    sigmahat2 = sum((y1 - lm(y1~xm21)$fitted.values)^2)/(n -half - length(ind1)-1)
  }
  sigmahat = (sigmahat1 + sigmahat2)/2
  return(list(sigma2 = sigmahat))
}
