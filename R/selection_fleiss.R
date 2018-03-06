# Fleiss: caclulate fleiss value for each beta corresponding to the lambda path and select best estimation
# input: complete_dataset(matrix): dataset used to estimate beta, first column is the value of y
#        lambda(vector): all the beta estimation needed to caculate BIC, each row is a set of
#                             estimation, there is no intercept term
#        repeat_b(scalar): time to repeat to get the average kappa
#        alpha_n(scalar): threshold to choose lambda, see definition in the paper
# output(list): kap_vec(vector): vector of all kappa value
#               VSS_min: value of lambda choosed by VSS
#


fleiss <- function(complete_data,cv.fold=NULL, cv.ind=NULL,alpha_n=0.1,penalty,lambda,gamma){
  #step 1 generate beta matrix
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
  flessi_vec <- flessi_kappa_calculator(lambda=lambda,logistic_list=logistic_list, penalty=penalty, gamma = gamma)

  #step2
  lambda_larger_idx <- which(flessi_vec/max(flessi_vec) >= 1-alpha_n)
  lambda_larger_value <- lambda[lambda_larger_idx]
  fleiss_min <- min(lambda_larger_value)


  # result <- list(selection_path = kappa_vec,
  #                selection_lambda_value = VSS_min)
  #
  result <- list(selection_path = flessi_vec,
                 selection_lambda_value = fleiss_min,
                 cv.ind=cv.ind)
  return(result)
}


############################################
#fleiss_kappa_calculator: calculate kappa value based on the input lambda(vector) for k fold
#input: lambda, penalty
#ouput: matrix with lambda and BIC value
flessi_kappa_calculator <- function(lambda,logistic_list, penalty, gamma){

  beta_list <- lapply(logistic_list,beta_est_path,penalty=penalty,lambda=lambda,gamma=gamma)
  length_min <- min(sapply(beta_list,nrow))
  beta_kappa_list <- list() #each element in this list is an matrix of beta used to caculate flessi kappa,there is NO intercept
  for(i in 1:length_min){
    beta_kappa_list[[i]] <- t(sapply(beta_list,FUN = function(matrix1) matrix1[i,]))[,-1]
  }
  kappa_vec <- sapply(beta_kappa_list,fleiss_kappa)
  return(kappa_vec)
}


# fleiss_kappa: define a function to calculate the fleiss for k beta_estimation
# input: beta_matrix(matrix): each row is a beta estimation, beta here should without intercept
# output(scalar): the kappa value

fleiss_kappa<- function(beta_matrix){
  binary_matrix <- apply(beta_matrix,2,FUN=function(row) as.numeric(row!=0))
  count_selected <- colSums(binary_matrix)
  confusion_matirx <- cbind(count_selected,nrow(beta_matrix)-count_selected)
  colnames(confusion_matirx) <- c("selection_count","unselection_count")
  if(sum(confusion_matirx[,1])==0 | sum(confusion_matirx[,2])==0){
    kappa <- -1
  }else{
    pe_vec <- apply(confusion_matirx, 2, sum)
    pe_vec <- pe_vec/(nrow(confusion_matirx)*nrow(beta_matrix))
    pe <- sum(pe_vec^2)
    p_vec <- apply(confusion_matirx, 1, FUN=function(x) sum(x*(x-1)))
    p_vec <- p_vec/(nrow(beta_matrix)*(nrow(beta_matrix)-1))
    p <- mean(p_vec)

    kappa <- (p-pe)/(1-pe)
  }
  return(kappa)
}
