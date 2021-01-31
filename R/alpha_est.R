#' Estimate alpha
#'
#' Estimate the value of best alpha that will be used for calculation
#' @param X a data frame
#' @param Y a vector of response variable
#' @param XXX_list a
#' @param rsquare
#' @param samp_time
#' @param n_var
#'
#' @return
#' @export
#'
#' @examples
alpha_estimate <- function(X, Y, XXX_list,rsquare,samp_time,n_var=0){
  if(n_var==0)
  {
    X1 <- scale(X)
    lambda_min <- gcv_lambda_min(X,Y)
    n_var <- sum(abs(glmnet(X, Y, lambda=lambda_min)$beta)>0)
  }
  # Estimate
  alpha_est <- 0
  for (i in 1:samp_time)
  {
    print(i)
    inf_index <- sample(1:ncol(X), n_var)
    Y <- Y_sim(X, inf_index, rsquare,mode="uniform")
    pr_auc <- 0
    for(j in 1:length(XXX_list))
    {
      Ecars_score <- abs(XXX_list[[j]] %*% scale(Y))/n
      pr_auc[j] <- comp_auc(as.vector(Ecars_score),inf_index,p=1)
    }
    alpha_est[i] <- alpha[which.max(pr_auc)]
  }
  return(alpha_est)
} #estimate alpha in ECAR
