#' fit and select variable(s) for data with missing value
#'
#' Fit a model based on a pseudo likelihood and select variable(s) through one of multiple techniques. The regularization path is computed for lasso, SCAD, or MCP.
#' Three steps are used to finish this the variable selection purpose: 1. remove missing and pair each observations;
#' 2. use penalty to get lambda path and corresponding beta matrix; 3. use specific method to finish variable selection.
#' @param x the covariate matrix, should be in matrix format and at least two columns, each row is an observation
#' @param y the response variable
#' @param penalty the penalty used for regularization, can be lasso, SCAD, or MCP. The default is lasso.
#' @param method the variable selection method, can be cross-validation (CV), Bayesian information criterion (BIC),
#' BIC1 and BIC2 are adapted for the consistency in the high dimension, sBIC is the information stability,
#' sBIC1 and sBIC2 are information stability for high dimension data, sVS is the variable selection stability,
#' sEST is the estimation stability
#' @param lambda lambda path used in the regularization path. If not specified by user, the path will be generated automatically
#' @param fold the number of folds used to divided data, will be used in CV, sBIC, sBIC1, sBIC2, sVS, and sEST method
#' @param cv.ind a vector to indicate what fold each observations belong, useful to make reproducible research
#' @param repeat_b B parameter in sVS method, the repeating time to calculate selection stability criteria
#' @param alpha_n the parameter used to take care of variables with weak effect in sVS method
#' @param refit If TRUE, refit technique will be used to get estimation, i.e., use selection variable to refit
#' the model to get estimation
#' @param use.penalty If TRUE, use penalty and variable selection techniques; if FALSE, just fit a logistic regression model with
#' paired data
#' @param gamma the tuning parameter of the SCAD/MCP. Default is 3.7 for SCAD and 3 for MCP
#'
#' @examples
#' n <- 50
#' p <- 8
#' beta <- c(3,0,1.5,0,2,rep(0,p-5))
#' xm <- matrix(rnorm(n*p),ncol = p, nrow = n)
#' y <- xm %*% beta + rnorm(n)
#' colnames(xm) <- paste0("Var_",1:p)
#'
#' fit01 <- tvsmiss(x=xm,y=y)
#' fit01$selection_beta
#'
#' fit02 <- tvsmiss(x=xm,y=y,method = "BIC")
#' fit02$selection_beta
#' fit02$beta_matrix
#'
#' fit06 <- tvsmiss(x=xm,y=y,penalty = "SCAD",method = "sVS",fold = 5)
#' fit06$selection_beta
#' fit06$beta_matrix
#'
#' @import glmnet
#' @importFrom stats approx binomial coef complete.cases glm model.matrix predict residuals sd
#' @importFrom utils combn
#' @useDynLib TVsMiss, .registration=TRUE
#'
#' @export
tvsmiss <- function(x,y,penalty=c("lasso", "MCP", "SCAD"),
                    method=c("CV", "BIC", "BIC1", "BIC2", "sBIC", "sBIC1", "sBIC2","sVS", "sEST"),
                    lambda=NULL,fold=5,cv.ind=NULL,repeat_b=20,alpha_n=0.1,refit=F,
                    gamma=switch(penalty, SCAD=3.7, MCP=3,lasso=NA),use.penalty=T){
  # three steps:
  # 1. remove missing and pair each observations
  # 2. use penalty to get lambda path and corresponding beta matrix
  # 3. use specific method to finish variable selection
  #
  # If cv.null is not NULL, then fold will be ignored
  penalty <- match.arg(penalty)
  method <- match.arg(method)
  this.call = match.call()

  if (class(x) != "matrix") {
    tmp <- try(x <- model.matrix(~0+., data=x), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("x must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(x)=="integer") storage.mode(x) <- "double"
  if (class(y) != "numeric") {
    tmp <- try(y <- as.numeric(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must numeric or able to be coerced to numeric")
  }

  if (length(y) != nrow(x)) stop("x and y do not have the same number of observations")

  #step 1: remove missing
  colnames(x) <- if (is.null(colnames(x))) paste("V",1:ncol(x),sep="") else colnames(x)

  data <- cbind(y,x)
  complete_idx <- c(1:length(y))[complete.cases(data)]
  y_complete <- y[complete_idx]
  x_complete <- x[complete_idx,]
  #generate logistic data
  sample_no_missing <- cbind(y_complete,x_complete)

  logistic_sample <- pairdata(sample_no_missing)

  if(use.penalty){

    if(is.null(cv.ind)){
      # check fold is in the correct range
      if(fold<2 | fold>nrow(sample_no_missing)) stop("fold should be greater than 1 and less than the rows in the complete data(after deleting missing)")
      cv.ind <- ceiling(sample(1:nrow(sample_no_missing))/(nrow(sample_no_missing)+sqrt(.Machine$double.eps))*fold)
    }else{
      if(length(cv.ind) != nrow(sample_no_missing) | max(cv.ind) > nrow(sample_no_missing)) stop("cv.ind is not match to the complete data")
    }

    sample_no_missing_list <- list()
    for(i in 1:max(cv.ind)){
      sample_no_missing_list[[i]] <- sample_no_missing[cv.ind == i,]
    }

    # logistic_list_estimation use (k-1) fold data to generate logistic dataset and estimate parameter
    # logistic_list_verification use 1 fold data to genete logistic dataset and calculate log-likelihood
    #if (method %in% c("CV", ))
    logistic_list_estimation <- cv_logistic_prepare(complete_data=sample_no_missing, cv.ind=cv.ind)
    # logistic_list_verification <- lapply(sample_no_missing_list, pairdata)

    #step2 use penatly to estimate
    current_model <- model_est_path(logistic_sample=logistic_sample,lambda=lambda,penalty=penalty,gamma=gamma)
    beta_matrix <- t(as.matrix(coef(current_model)))[,-1] #  remove intercept 0
    lambda <- current_model$lambda

    #step3 variable selections
    if(method == "CV"){
      begin_time <- proc.time()

      selection_res <- cv_sel(logistic_list=logistic_list_estimation, complete_data_list=sample_no_missing_list,
                              cv.fold=fold, cv.ind=cv.ind, penalty=penalty, lambda=lambda, gamma=gamma)
      selection_idx <- selection_res$lambda_idx
      selection_lambda <- lambda[selection_res$lambda_idx]
      selection_path <- selection_res$selection_path
      selection_cv.ind <- selection_res$cv.ind
      selection_beta <- beta_matrix[selection_idx,]
      final_time <- proc.time() - begin_time

    }else if(method == "BIC"){
      begin_time <- proc.time()

      selection_res <- BIC_log_cal(complete_dataset=sample_no_missing,beta_matrix=beta_matrix)
      selection_idx <- selection_res$lambda_idx
      selection_lambda <- lambda[selection_res$lambda_idx]
      selection_path <- selection_res$selection_path
      selection_cv.ind <- NULL
      selection_beta <- beta_matrix[selection_idx,]

      final_time <- proc.time() - begin_time

    }else if(method == "sBIC"){
      begin_time <- proc.time()

      selection_res <- sBIC(complete_data=sample_no_missing, cv.fold=fold, cv.ind=cv.ind,
                             penalty=penalty, lambda=lambda, gamma=gamma)
      selection_idx <- selection_res$lambda_idx
      selection_lambda <- lambda[selection_res$lambda_idx]
      selection_path <- selection_res$selection_path
      selection_cv.ind <- selection_res$cv.ind
      selection_beta <- beta_matrix[selection_idx,]

      final_time <- proc.time() - begin_time

    }else if(method == "sBIC1"){
      begin_time <- proc.time()

      selection_res <- sBIChigh(complete_data=sample_no_missing, cv.fold=fold, cv.ind=cv.ind,
                            penalty=penalty, lambda=lambda, gamma=gamma)
      selection_idx <- selection_res$lambda_idx
      selection_lambda <- lambda[selection_res$lambda_idx]
      selection_path <- selection_res$selection_path
      selection_cv.ind <- selection_res$cv.ind
      selection_beta <- beta_matrix[selection_idx,]

      final_time <- proc.time() - begin_time

    }else if(method == "sBIC2"){
      begin_time <- proc.time()

      selection_res <- sBICultrahigh(complete_data=sample_no_missing, cv.fold=fold, cv.ind=cv.ind,
                            penalty=penalty, lambda=lambda, gamma=gamma)
      selection_idx <- selection_res$lambda_idx
      selection_lambda <- lambda[selection_res$lambda_idx]
      selection_path <- selection_res$selection_path
      selection_cv.ind <- selection_res$cv.ind
      selection_beta <- beta_matrix[selection_idx,]

      final_time <- proc.time() - begin_time

    }else if(method == "sVS"){
      if(fold == 2){
        begin_time <- proc.time()

        selection_res <- VSS(complete_data=sample_no_missing,lambda=lambda,repeat_b=repeat_b,
                             alpha_n=alpha_n,penalty=penalty,gamma=gamma)
        selection_idx <- which(lambda == selection_res$selection_lambda_value)
        selection_lambda <- selection_res$selection_lambda_value
        selection_path <- selection_res$selection_path
        selection_cv.ind <- NULL
        selection_beta <- beta_matrix[selection_idx,]

        final_time <- proc.time() - begin_time
      }else{
        begin_time <- proc.time()

        selection_res <- fleiss(complete_data=sample_no_missing, cv.fold=fold, cv.ind=cv.ind,
                                penalty=penalty, lambda=lambda, alpha_n=alpha_n,gamma=gamma)
        selection_idx <- which(lambda == selection_res$selection_lambda_value)
        selection_lambda <- selection_res$selection_lambda_value
        selection_path <- selection_res$selection_path
        selection_cv.ind <- selection_res$cv.ind
        selection_beta <- beta_matrix[selection_idx,]

        final_time <- proc.time() - begin_time
      }
    }else if(method == "sEST"){
      begin_time <- proc.time()

      selection_res <- sEST(complete_data=sample_no_missing, cv.fold=fold, cv.ind=cv.ind,
                            penalty=penalty, lambda=lambda, gamma=gamma)
      selection_idx <- selection_res$lambda_idx
      selection_lambda <- lambda[selection_res$lambda_idx]
      selection_path <- selection_res$selection_path
      selection_cv.ind <- selection_res$cv.ind
      selection_beta <- beta_matrix[selection_idx,]

      final_time <- proc.time() - begin_time

    }else if(method == "BIC1"){
      begin_time <- proc.time()

      selection_res <- BIC_high(complete_dataset=sample_no_missing,beta_matrix=beta_matrix)
      selection_idx <- selection_res$lambda_idx
      selection_lambda <- lambda[selection_res$lambda_idx]
      selection_path <- selection_res$selection_path
      selection_cv.ind <- NULL
      selection_beta <- beta_matrix[selection_idx,]

      final_time <- proc.time() - begin_time
    }else if(method == "BIC2"){
      begin_time <- proc.time()

      selection_res <- BIC_ultrahigh(complete_dataset=sample_no_missing,beta_matrix=beta_matrix)
      selection_idx <- selection_res$lambda_idx
      selection_lambda <- lambda[selection_res$lambda_idx]
      selection_path <- selection_res$selection_path
      selection_cv.ind <- NULL
      selection_beta <- beta_matrix[selection_idx,]

      final_time <- proc.time() - begin_time
    }

  }else{
    begin_time <- proc.time()

    glm_model <- glm(logistic_sample[,1]~-1+logistic_sample[,-1], family = binomial(link = "logit"))
    current_model = glm_model
    beta_matrix = NULL
    lambda = NULL
    cv.ind = NULL
    selection_idx = NULL
    selection_lambda = NULL
    selection_path = NULL
    selection_cv.ind = NULL
    selection_beta =  current_model$coefficients
    names(selection_beta) <- colnames(logistic_sample)[-1]

    final_time <- proc.time() - begin_time
  }

  deviance.ratio = 1 + loglikelihood_r(c(0,selection_beta),sample_no_missing)/log(2)

  #refit
  if(refit & use.penalty){
    if(sum(selection_beta != 0) == 0){
      refit_beta = selection_beta
    }else{
      refit_model <- glm(logistic_sample[,1]~-1+logistic_sample[,-1][,which(selection_beta != 0)], family = binomial(link = "logit"))
      refit_beta <- rep(0,ncol(logistic_sample)-1)
      refit_beta[which(selection_beta != 0)] <- refit_model$coefficients
      names(refit_beta) <- colnames(logistic_sample)[-1]
    }


  }else{
    refit_beta=NULL
  }


  res <- list(ls=logistic_sample,c_idx=complete_idx,model=current_model,
              beta_matrix=beta_matrix,lambda=lambda, cv.ind=cv.ind, fold=fold,
              selection_idx = selection_idx, selection_lambda = selection_lambda,
              selection_path = selection_path, selection_cv.ind = selection_cv.ind,
              selection_beta =  selection_beta,
              refit_beta = refit_beta,
              null.deviance = 2*log(2),
              deviance.ratio=deviance.ratio,
              gamma=gamma,
              call = this.call,
              running_time = final_time)
  class(res) = "TVsMiss"
  return(res)
}
