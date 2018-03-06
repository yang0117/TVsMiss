############################################
#cv_logistic_prepare: generate the logistic dataset in a list for cv, use (k-1) fold of data each time
#the purpose of this function is save time by avoiding generate logistic list every time
#input: sample_original,cv.ind
#output: a list with all logistic dataset for cross valiadation(our method)

cv_logistic_prepare <- function(complete_data,cv.ind){
  #generate the logistic dataset in a list
  logistic_list <- list()
  for (i in 1:max(cv.ind)){
    data_k <- complete_data[cv.ind!=i,]
    logistic_list[[i]] <- pairdata(data_k)
  }
  return(logistic_list)
}

