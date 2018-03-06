# sBIC: caculate BIC value with "BIC stability", it is like CV method, use (k-1)/k portion
#        data to calculate and use one protion data to get BIC value and then get the mean for k BIC value

# input: original_dataset(matrix): dataset used to estimate beta, first column is the value of y
#        lambda_path(vector): all the beta estimation needed to caculate BIC, each row is a set of
#                             estimation, there is no intercept term
#        k_BIC: the fold used for ESCV in cross validation
#         logistic_list(list): a list of logistic data list based on K-fold prepared for BIC to save time
#
# output(list): BIC_stable_est(vector): vector of all BIC value
#               BIC_stable_value: value of lambda choosed by BIC_stable method

sBIChigh<- function(complete_data, cv.fold=NULL, cv.ind=NULL, penalty, lambda, gamma){
  # check cv.fold is in the correct range
  if(cv.fold<2 | cv.fold>nrow(complete_data)) stop("cv.fold should be greater than 1 and less than the rows in the complete data(after deleting missing)")

  if(is.null(cv.ind)){
    cv.ind <- ceiling(sample(1:nrow(complete_data))/(nrow(complete_data)+sqrt(.Machine$double.eps))*cv.fold)
  }

  complete_data_list <- list()
  for(i in 1:max(cv.ind)){
    complete_data_list[[i]] <- complete_data[cv.ind == i,]
  }
  logistic_list <- cv_logistic_prepare(complete_data=complete_data,cv.ind=cv.ind)

  #create a matrix to store lambda value and corresponding BIC value
  path_matrix <- sBIChigh_calculator(lambda = lambda, complete_data = complete_data,
                                logistic_list = logistic_list, complete_data_list = complete_data_list,
                                penalty = penalty, gamma=gamma)

  result <- list(lambda_idx = which.min(path_matrix[,2]),
                 selection_path = path_matrix[,2],
                 cv.ind=cv.ind)
  return(result)
}



############################################
#BIC_calculator: calculate BIC value based on the input lambda(vector) for k fold
#input: lambda, complete_data, penalty
#ouput: matrix with lambda and BIC value
sBIChigh_calculator <- function(lambda, complete_data, logistic_list, complete_data_list, penalty,gamma){

  beta_list <- lapply(logistic_list,beta_est_path,penalty=penalty,lambda=lambda,gamma=gamma)
  res_list <- list()

  for(i in 1:length(beta_list)){
    res_list[[i]] <- apply(beta_list[[i]], 1, BIC_high_singel_value, dataset=complete_data_list[[i]])
  }

  length_min <- min(sapply(res_list,length))
  cv_res_matrix <- matrix(-999,length_min,length(beta_list))
  for(i in 1:length(res_list)){
    cv_res_matrix[,i] <- res_list[[i]][1:length_min]
  }
  result <- cbind(lambda[1:length_min],rowMeans(cv_res_matrix))
  colnames(result) <- c(paste("lambda_",penalty,sep = ""),"BIC")

  return(result)
}
