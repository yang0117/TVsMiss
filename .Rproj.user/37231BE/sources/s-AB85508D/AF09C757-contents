########################################################
#function to calculate the log-pseudolikelihood
#input:dataset, beta_estimation HAVE intercept
#output: likelihood value
loglikelihood_r<- function(beta_estimation,dataset){
  #define a function to calculate each row of new logistic dataset
  single_log_likelihood <- function(idx_vec,dataset_y,dataset_x,beta_estimation){
    res <- (-log(1+exp(-(dataset_y[idx_vec[1]]-dataset_y[idx_vec[2]])*
                         sum(as.vector(beta_estimation)*as.vector(dataset_x[idx_vec[1],]-dataset_x[idx_vec[2],])))))
    return(res)
  }

  #prepare data
  dataset <- as.matrix(dataset)
  row_idx <- combn(dim(dataset)[1],2)
  dataset_y <- dataset[,1]
  #wihtout 1 when there is no intercept
  dataset_x <- cbind(1,dataset[,-1])
  # print(as.vector(dataset_x[1,]-dataset_x[1,]))
  # print(as.vector(beta_estimation))
  if(length(as.vector(dataset_x[1,]-dataset_x[1,])) != length(as.vector(beta_estimation))) stop("x and beta_est does not match loglikelihood")

  # print("begin")
  # time1 <- proc.time()
  loglikelihood_vec <- apply(row_idx,2,single_log_likelihood,dataset_y=dataset_y,dataset_x=dataset_x,beta_estimation=beta_estimation)
  # time1 <- proc.time() - time1
  # print(time1)
  # print("done")

  loglikelihood = mean(loglikelihood_vec)

  return(loglikelihood)
}


# ########################################################
# #function to calculate the log-pseudolikelihood for logistic data
# #input:dataset, beta_estimation HAVE intercept
# #output: likelihood value
# loglikelihood_logistic <- function(beta_estimation,dataset){
#
#   single_log_likelihood <- function(data_one_obs, beta_estimation){
#     x <- as.vector(c(1,data_one_obs[-1]))
#     x_beta <- sum(x*as.vector(beta_estimation))
#     res <- log(exp(data_one_obs[1]*x_beta)/(1+exp(x_beta)))
#     return(res)
#   }
#
#   # print(dim(dataset))
#   # print("begin")
#   # time1 <- proc.time()
#   loglikelihood_vec <- apply(dataset,1,single_log_likelihood,beta_estimation=beta_estimation)
#   # time1 <- proc.time() - time1
#   # print(time1)
#   # print("done")
#
#   loglikelihood = mean(loglikelihood_vec)
#
#   return(loglikelihood)
# }




