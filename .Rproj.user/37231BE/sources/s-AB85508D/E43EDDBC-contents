# sEST: caclulate sEST value for each beta corresponding to the lambda path and select best estimation
# input: complete_data(matrix): dataset used to estimate beta, first column is the value of y
#        lambda_path(vector): all the beta estimation needed to caculate BIC, each row is a set of
#                             estimation, there is no intercept term
#        cv.fold: the fold used for sEST in cross validation
#        lambda_value_ll: lambda value from the cv based on log-likelihood
# output(list): sEST_mse_lambda: lambda value chosen by cv MSE method,
#               sEST_ll_lambda: lambda value chosen by cv ll method,
#               ES_vec: ES value vector

sEST <- function(complete_data, lambda, cv.fold=NULL, cv.ind=NULL, penalty,gamma){
  # step1: divide data into five parts
  # step2: get estimation for each beta
  # step3: calculate sEST
  # step4: choose lambda position

  # step1
  # generate a list to store the obs number be removed in cross validation
  # check cv.fold is in the correct range
  if(cv.fold<2 | cv.fold>nrow(complete_data)) stop("cv.fold should be greater than 1 and less than the rows in the complete data(after deleting missing)")

  if(is.null(cv.ind)){
    cv.ind <- ceiling(sample(1:nrow(complete_data))/(nrow(complete_data)+sqrt(.Machine$double.eps))*cv.fold)
  }

  #generate the original and logistic dataset for cross validation, each dataset contains (k-1)/k data
  complete_data_list <- list()
  for(i in 1:max(cv.ind)){
    complete_data_list[[i]] <- complete_data[cv.ind == i,]
  }
  logistic_list <- cv_logistic_prepare(complete_data=complete_data,cv.ind=cv.ind)
  # step2
  # this beta HAS intercept
  beta_list <- lapply(logistic_list,beta_est_path,penalty=penalty,lambda=lambda,gamma=gamma)
  beta_min_length <- min(sapply(beta_list,nrow))
  beta_list <- lapply(beta_list,FUN=function(x) x[1:beta_min_length,])

  #Note: in the following steps, y_hat is based on complete data

  #step3
  y_hat_list <- list()
  logistic_data <- pairdata(complete_data)
  for(i in 1:length(beta_list)){
    y_hat_list[[i]] <- as.matrix(logistic_data[,-1]) %*% t(beta_list[[i]][,-1])
  }

  y_hat_full_matrix <- Reduce('+',y_hat_list)/length(y_hat_list)
  l2_diff_matrix <- matrix(0,nrow=max(cv.ind),ncol = beta_min_length)
  for(i in 1:length(y_hat_list)){
    l2_diff_matrix[i,] <- apply((y_hat_list[[i]] - y_hat_full_matrix)^2, 2, sum)
  }
  l2_diff_mean <- apply(l2_diff_matrix,2,mean)
  ES_vec <- l2_diff_mean/apply(y_hat_full_matrix^2,2,sum)

  #step4
  ES_grad = ES_vec[2:length(ES_vec)]- ES_vec[1:(length(ES_vec)-1)]
  ES_neg_idx = which(ES_grad < 0)
  if(length(ES_neg_idx)==0) {
    ES_start_idx = 1 #this value will be shared by both MSE and log-likelihood method
  } else {ES_start_idx = ES_neg_idx[1]}

  #get lambda value from cv
  cv_res <- cv_sel(logistic_list=logistic_list, complete_data_list=complete_data_list,
                   cv.fold=cv.fold, cv.ind=cv.ind,
                   penalty=penalty, lambda=lambda,gamma=gamma)
  cv_idx <- cv_res$lambda_idx

  #find ES idx based on cv from log-likelihood
  if(ES_start_idx > cv_idx) {
    ES_ll_start_idx = cv_idx
  }else{ES_ll_start_idx=ES_start_idx}
  sEST_ll_idx = which.min(ES_vec[ES_ll_start_idx:cv_idx]) + ES_ll_start_idx - 1

  res <- list(lambda_idx = sEST_ll_idx,
              selection_path = ES_vec,
              cv.ind=cv.ind)
  return(res)
}
