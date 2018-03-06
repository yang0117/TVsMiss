#CV method

cv_sel <- function(logistic_list, complete_data_list, cv.fold=NULL, cv.ind=NULL, penalty, lambda, gamma){

  # # check cv.fold is in the correct range
  # if(cv.fold<2 | cv.fold>nrow(complete_data)) stop("cv.fold should be greater than 1 and less than the rows in the complete data(after deleting missing)")
  #
  # if(is.null(cv.ind)){
  #   cv.ind <- ceiling(sample(1:nrow(complete_data))/(nrow(complete_data)+sqrt(.Machine$double.eps))*cv.fold)
  # }
  #
  # complete_data_list <- list()
  # for(i in 1:max(cv.ind)){
  #   complete_data_list[[i]] <- complete_data[cv.ind == i,]
  # }
  # logistic_list <- cv_logistic_prepare(complete_data=complete_data,cv.ind=cv.ind)

  #create a matrix to store lambda value and corresponding cv value
  path_matrix <- cv_calculator(lambda = lambda, logistic_list = logistic_list,
                               complete_data_list = complete_data_list, penalty = penalty, gamma=gamma)

  result <- list(lambda_idx = which.max(path_matrix[,2]),
                 selection_path = path_matrix[,2],
                 cv.ind=cv.ind)
  return(result)
}

############################################
#cv_calculator: calculate cv values(our method) based on the input lambda(scalar)
#input: lambda, sample_original, k, penalty_indicator
#ouput: a single cv_value corresponding to lambda and penalty_indicator
cv_calculator <- function(lambda,logistic_list,complete_data_list,penalty,gamma){

  cv_beta_list <- list()

  time1 <- proc.time()
  # for(i in 1:length(logistic_list)){
  #   print(dim(logistic_list[[i]]))
  #   cv_beta_list[[i]] <- beta_est_path(logistic_list[[i]],penalty=penalty,lambda=lambda)
  # }
  cv_beta_list <- lapply(logistic_list,beta_est_path,penalty=penalty,lambda=lambda,gamma=gamma)
  time1 <- proc.time() - time1
  #print("beta estimation time")
  #print(time1)

  cv_res_list <- list()

  loglikelihood_c <- function(beta,complete_data){
    res <- .Call("loglikelihood",beta,complete_data)
    return(res)
  }

  time1 <- proc.time()
  for(i in 1:length(cv_beta_list)){
    #cv_res_list[[i]] <- apply(cv_beta_list[[i]],1,loglikelihood_logistic,complete_data_list[[i]])
    cv_res_list[[i]] <- apply(cv_beta_list[[i]],1,loglikelihood_c,complete_data_list[[i]])
  }
  time1 <- proc.time() - time1
  #print("log-likelihood time")
  #print(time1)

  length_min <- min(sapply(cv_res_list,length))
  cv_res_matrix <- matrix(-999,length_min,length(cv_beta_list))
  for(i in 1:length(cv_res_list)){
    cv_res_matrix[,i] <- cv_res_list[[i]][1:length_min]
  }
  result <- cbind(lambda[1:length_min],rowMeans(cv_res_matrix))
  colnames(result) <- c(paste("lambda_",penalty,sep = ""),"cv_value")
return(result)
}
