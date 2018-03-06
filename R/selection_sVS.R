# VSS: caclulate VSS value for each beta corresponding to the lambda path and select best estimation
# input: complete_dataset(matrix): dataset used to estimate beta, first column is the value of y
#        lambda(vector): all the beta estimation needed to caculate BIC, each row is a set of
#                             estimation, there is no intercept term
#        repeat_b(scalar): time to repeat to get the average kappa
#        alpha_n(scalar): threshold to choose lambda, see definition in the paper
# output(list): kap_vec(vector): vector of all kappa value
#               VSS_min: value of lambda choosed by VSS


VSS <- function(complete_data,lambda,repeat_b=20,alpha_n=0.1,penalty,gamma){
  # function work as step 1 and step 2 in the paper, current_lambda is a single value

  #step1 split data and estimate
  split_idx <- t(replicate(repeat_b,ceiling(sample(1:nrow(complete_data))/(nrow(complete_data)+sqrt(.Machine$double.eps))*2)))
  #a list of size repeat b
  data_splitted <- apply(split_idx,1,FUN = function(s_idx,data) split.data.frame(data,s_idx), data=complete_data)
  # data_a <- lapply(data_splitted, FUN=function(x) x[[1]])
  # data_b <- lapply(data_splitted, FUN=function(x) x[[2]])
  data_a_logistic <- lapply(data_splitted, FUN=function(x) pairdata(x[[1]]))
  data_b_logistic <- lapply(data_splitted, FUN=function(x) pairdata(x[[2]]))

  data_a_est <- lapply(data_a_logistic, beta_est_path, lambda=lambda, penalty=penalty, gamma=gamma)
  data_b_est <- lapply(data_b_logistic, beta_est_path, lambda=lambda, penalty=penalty, gamma=gamma)

  #step2/step3 caculate kappa
  #remove intercept
  kappa_res <- matrix(-99,repeat_b,length(lambda))

  for(i in 1:repeat_b){
    current_a <- data_a_est[[i]]
    current_b <- data_b_est[[i]]
    if(nrow(current_a) > nrow(current_b)){
      current_b = rbind(current_b, t(replicate(nrow(current_a) - nrow(current_b), current_b[nrow(current_b),])))
    }else if(nrow(current_a) < nrow(current_b)){
      current_a = rbind(current_a, t(replicate(nrow(current_b) - nrow(current_a), current_a[nrow(current_a),])))
    }

    current_kappa_vec <- sapply(1:nrow(current_a), FUN=function(row_idx) kappa_cal(current_a[row_idx,-1], current_b[row_idx,-1]))
    if(length(current_kappa_vec) < length(kappa_res[i,])){
      tail_vec <- rep(current_kappa_vec[length(current_kappa_vec)],length(kappa_res[i,])-length(current_kappa_vec))
      current_kappa_vec <- c(current_kappa_vec,tail_vec)
    }
    kappa_res[i,] <- current_kappa_vec

    # kappa_res[i,] <- sapply(1:nrow(data_a_est[[i]]), FUN=function(row_idx,list_idx) kappa_cal(data_a_est[[list_idx]][row_idx,-1], data_b_est[[list_idx]][row_idx,-1]),
    #                         list_idx=i)
  }

  kappa_vec <- colMeans(kappa_res)

  #step4
  lambda_larger_idx <- which(kappa_vec/max(kappa_vec) >= 1-alpha_n)
  if(length(lambda_larger_idx) == 0) stop('none of lambda match the threshold(07.03.VSS/VSS)')
  lambda_larger_value <- lambda[lambda_larger_idx]
  VSS_min <- min(lambda_larger_value)


  result <- list(selection_path = kappa_vec,
                 selection_lambda_value = VSS_min)
  return(result)
}


# kappa_cal: define a function to calculate the kappa defined in the variable selection stability paper
# input: beta_est1, beta_est2: two estimation from splitted dataset,beta here should without intercept
# output(scalar): the kappa value

kappa_cal <- function(beta_est1,beta_est2){
  if (length(beta_est1) != length(beta_est2)) stop("lenght of two beta does not match(07.03.VSS/kappa_cal)")
  #get non-zero postion
  beta_est1_nonzero <- which(beta_est1 != 0)
  beta_est2_nonzero <- which(beta_est2 != 0)
  beta_est1_nonzero_c <- which(beta_est1 == 0)
  beta_est2_nonzero_c <- which(beta_est2 == 0)
  p = length(beta_est1)

  #calculate
  n11 = length(intersect(beta_est1_nonzero,beta_est2_nonzero))
  n12 = length(intersect(beta_est1_nonzero,beta_est2_nonzero_c))
  n21 = length(intersect(beta_est1_nonzero_c,beta_est2_nonzero))
  n22 = length(intersect(beta_est1_nonzero_c,beta_est2_nonzero_c))

  if(length(beta_est2_nonzero_c) == length(beta_est2_nonzero_c) & length(beta_est2_nonzero_c) == p){
    kap <- -1
  }else if(length(beta_est2_nonzero) == length(beta_est2_nonzero) & length(beta_est2_nonzero) == p){
    kap <- -1
  }else if(n12 == n21 & n12 == 0){
    kap <- 1
  }else if(n11 == n22 & n11 == 0 & n12 == n21 & n12 == p/2){
    kap <- -1
  }else{
    pra = (n11 + n22)/p
    pre = (n11 + n12)*(n11 + n21)/p^2 + (n12 + n22)*(n21 + n22)/p^2

    kap <- (pra - pre)/(1 - pre)
  }
  return(kap)
}
