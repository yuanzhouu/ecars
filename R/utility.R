
# load functions

gcv_lambda_min <- function(X,Y,return="lambda"){
  ls.fit <- glmnet(X,Y,family="gaussian")
  df <- ls.fit$df
  df <- df[df<nrow(X)]
  res <- data.frame(df=df)
  for(i in 1:nrow(res))
  {
    res$SSE[i] <- sum((Y-ls.fit$a0[i]-X %*% ls.fit$beta[,i])^2)
    res$gcv_SSE[i] <- res$SSE[i]/nrow(X)/(1-res$df[i]/nrow(X))^2
  }
  if(return=="lambda")
  {return(ls.fit$lambda[which.min(res$gcv_SSE)])
  }else if(return=="index")
  {
    return(which.min(res$gcv_SSE))
  }
} # gcv-lasso
gcv_lambda_ridge_min <- function(X, Y){
  ridge.fit <- glmnet(X,Y,family="gaussian",alpha=0)
  v_lambda <- svd(X)$d^2
  v_lambda[is.na(v_lambda)] <- 0
  rg_res <- data.frame(lambda=ridge.fit$lambda)
  for(i in 1:100)
  {
    rg_res$ridge_trace[i] <- sum(v_lambda/(v_lambda+rg_res$lambda[i]*nrow(X)))
    rg_res$SSE[i] <- sum((scale(Y)-ridge.fit$a0[i]-X %*% ridge.fit$beta[,i])^2)
    rg_res$gcv_SSE[i] <- rg_res$SSE[i]/nrow(X)/(1-rg_res$ridge_trace[i]/nrow(X))^2
  }
  return(ridge.fit$lambda[which.min(res$gcv_SSE)])
} # gcv-ridge,under modification
replace_singular_X <- function(X){
  ind1 <- which(duplicated(t(X)))
  ind2 <- which(apply(X,2,var)==0)
  ind <- union(ind1,ind2)
  X[,ind] <- matrix(rnorm(length(ind)*nrow(X)),nrow(X),length(ind))
  return(X)
} # remove zero variance and duplicate X
get_stabls_set <- function(X, Y, q){
  res_set <- NULL
  for(k in 1:100)
  {
    #rand_id <- 1:nrow(X)
    rand_id <- sample(1:nrow(X), floor(nrow(X)/2))
    #lars.fit <- lars(X[rand_id,],Y[rand_id],type="lasso",max.steps = 2*floor(q),use.Gram = F)
    #sel_set <- which(!lars.fit$entry==0)
    ls.fit <- glmnet(X[rand_id,],Y[rand_id])
    if(all(ls.fit$df<q))
    {sel_set <- which(abs(ls.fit$beta[,min(ls.fit$lambda)])>0)
    }else
    {
      sel_set <- which(abs(ls.fit$beta[,which(ls.fit$df>q)[1]])>0)
    }
    res_set <- c(res_set,sel_set)
  }
  return(as.numeric(names(sort(table(res_set),decreasing = T))))
} # rank variables according to probability
get_sis_set <- function(X, Y){
  sis_score <- cor(X,Y)
  sis_score[is.na(sis_score)] <- 0
  return(sort.int(abs(sis_score),decreasing=T,index.return=T)$ix[1:ncol(X)])
} # rank variables acc to cor
get_ridge_set <- function(X,Y){


  ridge.fit_cv <- cv.glmnet(X,Y,alpha=0)
  ridge.fit <- glmnet(X,Y,alpha=0,lambda=ridge.fit_cv$lambda.min)
  return(sort.int(abs(as.vector(ridge.fit$beta)),
                  decreasing = T,index.return = T)$ix[1:ncol(X)])
} # rank varialbles acc to ridge coefficient
get_cars_set <- function(X, Y){
  corr <- cor(X,Y)
  corr[is.na(corr)] <- 0
  cars_score <- crossprod.powcor.shrink(X,corr,alpha=-0.5,verbose=F)
  return(sort.int(abs(cars_score),decreasing=T,index.return=T)$ix[1:ncol(X)])
} # rank varialbles acc to coefficient
get_Ecars_set <- function(X, Y, alpha){
  corr <- cor(X,Y)
  corr[is.na(corr)] <- 0
  Ecars_score <- crossprod.powcor.shrink(X,corr,alpha=alpha,verbose=F)
  return(sort.int(abs(Ecars_score),decreasing=T,index.return=T)$ix[1:ncol(X)])
} #rank varialbles acc to ridge coefficient
get_XXX_list <- function(X, alpha){
  X1 <- scale(X)
  XXX_list <- list()
  for(i in 1:length(alpha))
  {

    XXX_list[[i]] <- crossprod.powcor.shrink(X1,t(X1),alpha[i],verbose=F)
  }
  return(XXX_list)
} # (RX)^{-a}%*%X
Y_sim <- function(X,inf_index,rsquare,mode,abs=T){
  get_Y3 <- function(X, eta, ratio){
    sigma2 <- var(eta)*(1-ratio)/ratio
    lth <- nrow(X)
    Y <- rnorm(lth, eta, sd=sqrt(sigma2))
  }
  if(mode == "uniform")
  {
    coef <- runif(length(inf_index))
  }else if (mode == "normal")
  {
    coef <- rnorm(length(inf_index))
  }
  if(abs==T)
  {
    coef <- abs(coef)
  }
  eta <- apply(X[, inf_index], 1, function(x){
    return( sum(x*coef))
  }) #?պôﵽ25%???ͷ?????ˮƽ
  Y <- get_Y3(X, eta, ratio=rsquare)
  rm(.Random.seed, envir=.GlobalEnv)
  return(Y)
} # generate Y
comp_auc <- function(score, inf_index ,p){
  label_mat <- rep(0, length(score))
  label_mat[inf_index] <- 1
  method_res <- roc(label_mat,1-as.numeric(score),quiet=T)
  method_cords <- coords(method_res, "all", ret = c("recall", "precision"),
                         transpose=TRUE)
  lth <- ncol(method_cords)-1
  method_auc <- trapz(method_cords[1,lth:p],
                      method_cords[2,lth:p])
  return(method_auc)
} #compute PR-AUC
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
comp_pred_err <- function(X_train,Y_train,X_test,Y_test,sel_set,seq_num){
  err_seq <- 0
  for (i in 1:length(seq_num))
  {
    #????lasso
    ls.fit <- glmnet(X_train[,sel_set[1:seq_num[i]]],Y_train,family="gaussian",lambda=
                       gcv_lambda_min(X_train[,sel_set[1:seq_num[i]]],Y_train))
    pred <- predict(ls.fit,newx = X_test[,sel_set[1:seq_num[i]]])
    err_seq[i] <- sum((Y_test-pred)^2)/length(Y_test)
  }
  return(err_seq)
}
comp_pred_err_ridge <- function(X_train,Y_train,X_test,Y_test,sel_set,seq_num){
  err_seq <- 0
  for (i in 1:length(seq_num))
  {
    #????lasso
    ridge.fit_cv <- cv.glmnet(X_train[,sel_set[1:seq_num[i]]],Y_train,family="gaussian",
                              alpha=0)
    ridge.fit <- glmnet(X_train[,sel_set[1:seq_num[i]]],Y_train,family="gaussian",
                        alpha=0,lambda=ridge.fit_cv$lambda.min)
    pred <- predict(ridge.fit,newx = X_test[,sel_set[1:seq_num[i]]])
    err_seq[i] <- sum((Y_test-pred)^2)/length(Y_test)
  }
  return(err_seq)
}
VAR_RCV1 <- function (y, x) {
  n = nrow(x)
  p = ncol(x)
  half = ceiling(n/2)
  x1 = x[1:half, ]
  y1 = y[1:half]
  x2 = x[-c(1:half), ]
  y2 = y[-c(1:half)]
  fit.cv = cv.glmnet(x1, y1, family = "gaussian")
  ind = which.min(fit.cv$cvm)
  fits = fit.cv$glmnet.fit
  if (ind == 1) {
    sigmahat1 = sum(y2^2)/half
  }
  else {
    ind1 = which(abs(fits$beta[, ind]) > 0)
    if (length(ind1) > half/2) {
      betasort = sort(abs(fits$beta[ind1, ind]), decreasing = TRUE,
                      index.return = T)
      ind1 = ind1[betasort$ix[1:ceiling(half/2)]]
    }
    xm12 = x2[, ind1]
    sigmahat1 =  sum((y2 - lm(y2~xm12)$fitted.values)^2)/(half -
                                                            length(ind1))
  }
  fit.cv = cv.glmnet(x2, y2, family = "gaussian")
  ind = which.min(fit.cv$cvm)
  fits = fit.cv$glmnet.fit
  if (ind == 1) {
    sigmahat2 = sum(y1^2)/(n - half)
  }
  else {
    ind1 = which(abs(fits$beta[, ind]) > 0)
    if (length(ind1) > half/2) {
      betasort = sort(abs(fits$beta[ind1, ind]), decreasing = TRUE,
                      index.return = T)
      ind1 = ind1[betasort$ix[1:ceiling(half/2)]]
    }
    xm21 = x1[, ind1]

    sigmahat2 = sum((y1 - lm(y1~xm21)$fitted.values)^2)/(n -
                                                           half - length(ind1))
  }
  sigmahat = (sigmahat1 + sigmahat2)/2
  return(list(sigma2 = sigmahat))
} # estimate residual variance
get_cvlasso_score <- function(X, Y){
  ls_coef_cv <- glmnet(X, Y,alpha=1, lambda=gcv_lambda_min(X,Y))
  ls_cv_score <- rep(0,ncol(X))
  ls_cv_score[((ls_coef_cv$beta)@i)+1] <- abs((ls_coef_cv$beta)@x)
  ls_cv_score[-(((ls_coef_cv$beta)@i)+1)] <- runif(ncol(X)-length((ls_coef_cv$beta)@i),
                                                   0,min(abs((ls_coef_cv$beta)@x)))
  return(ls_cv_score)
} # compute score of selected, nonselected random
combine_list <- function(list1, list2){
  list3 <- list()
  for(i in 1:length(list1))
  {
    list3[[i]] <- rbind(list1[[i]], list2[[i]])
  }
  return(list3)
}
get_stabls_score <- function(X, Y, q){
  res_set <- NULL
  for(k in 1:100)
  {
    ls.fit <- glmnet(X,Y)
    if(all(ls.fit$df<q))
    {
      sel_set <- which(abs(ls.fit$beta[,min(ls.fit$lambda)])>0)
    }else
    {
      sel_set <- which(abs(ls.fit$beta[,which(ls.fit$df>q)[1]])>0)
    }
    res_set <- c(res_set,sel_set)
  }
  stabls_score <- rep(0,ncol(X))
  ind <- as.numeric(names(sort(table(res_set),decreasing = T)))
  stabls_score[ind] <- seq(1,0.5,length.out=length(ind))
  stabls_score[-ind] <- runif(ncol(X)-length(ind),
                              0,0.5)
  return(stabls_score)

} # compute stabls score
get_ridge_score <- function(X, Y){
  ridge.res <- cv.glmnet(X,Y,family="gaussian",alpha=0)
  ridge.res <- glmnet(X,Y,family="gaussian",alpha=0,
                      lambda=ridge.res$lambda.min)
  ridge_score <- abs((ridge.res$beta)@x)
  return(ridge_score)
} # compute ridge core
get_auc_meansd <- function(auc_mat,i,total){
  auc <- round(mean((auc_mat)[seq(i,length(auc_mat),by=total)]),digit=3)
  sd <- round(sd((auc_mat)[seq(i,length(auc_mat),by=total)]),digit=3)
  return(paste0(auc,"(",sd,")"))
}
comp_line <- function(lmscore,inf_index, n=300){
  lmscore_ind <- sort.int(lmscore,index.return = T,decreasing=T)$ix[1:n]
  lm_ind <- lmscore_ind %in% inf_index
  lm_ind <- cumsum(lm_ind)
  return(lm_ind)
} #compute the selection path
VAR_RCV2 <- function (y, x) {
  n = nrow(x)
  p = ncol(x)
  half = ceiling(n/2)
  x1 = x[1:half, ]
  y1 = y[1:half]
  x2 = x[-c(1:half), ]
  y2 = y[-c(1:half)]
  fit = glmnet(x1, y1, family = "gaussian")
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
  fit = glmnet(x2, y2, family = "gaussian")
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
