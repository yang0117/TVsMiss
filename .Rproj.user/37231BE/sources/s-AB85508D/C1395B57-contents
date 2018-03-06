# 07.02 BIC selection
# BIC_log_cal: this function used to caculate the BIC value of each group of beta estimation and
#          choose the estimation with minimimal BIC value, and the value log-likelihood is from logistics
# input: complete_dataset(matrix): dataset used to estimate beta, first column is the value of y
#        beta_matrix(matrix): all the beta estimation needed to caculate BIC, each row is a set of
#                             estimation, they NOT have intercept term
# output(list): BIC_est(vector): vector of all BIC value
#               BIC_idx: location of BIC choosed

BIC_log_cal <- function(complete_dataset,beta_matrix){
  if ((ncol(complete_dataset)-1) != (ncol(beta_matrix))) stop("number of beta does not math complete_dataset")
  BIC_vec <- apply(cbind(0,beta_matrix), 1, BIC_singel_value, dataset=complete_dataset)

  BIC_index <- which.min(BIC_vec)

  result <- list(lambda_idx = BIC_index,
                 selection_path = BIC_vec)
  return(result)
}

#define a function to cacluate the BIC value for single beta estimation
#single_beta_est HAVE the intercept
BIC_singel_value <- function(single_beta_est,dataset){
  # note, in the paper, there is no negative sign before likelihood function, this is because in
  # the paper l(beta) is the negative of log-likeloohd function
  # bic_value <- -2*loglikelihood(dataset=dataset,beta_estimation=single_beta_est) + sum(single_beta_est[-1] != 0)*log(choose(dim(dataset)[1],2))/choose(dim(dataset)[1],2)

  bic_value <- -2*.Call("loglikelihood",single_beta_est,dataset) + sum(single_beta_est[-1] != 0)*log(nrow(dataset))/nrow(dataset)
  return(bic_value)
}
