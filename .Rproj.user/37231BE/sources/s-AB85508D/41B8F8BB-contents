#input lambda path, penalty, return a model

model_est_path <- function(logistic_sample,lambda,penalty,gamma){
  if(penalty == "lasso"){
    current_model <- glmnet(as.matrix(logistic_sample[,-1]),factor(logistic_sample[,1],levels=c(0,1)),
                            family="binomial",intercept = F,standardize = F,lambda = lambda)
  }else if (penalty == "SCAD"){
    current_model <- suppressWarnings(ncvint(as.matrix(logistic_sample[,-1]),logistic_sample[,1],
                                             family="binomial",penalty="SCAD",lambda=lambda,gamma=gamma))
  }else if (penalty == "MCP"){
    current_model <- suppressWarnings(ncvint(as.matrix(logistic_sample[,-1]),logistic_sample[,1],
                                             family="binomial",penalty="MCP",lambda=lambda,gamma=gamma))
  }
  return(current_model)
}


#estimation result here HAVE intercept
#each row is a beta estimation
beta_est_path <- function(logistic_sample,lambda,penalty,gamma){
  if(penalty == "lasso"){
    current_model <- glmnet(as.matrix(logistic_sample[,-1]), factor(logistic_sample[,1],levels=c(0,1)),
                            family="binomial",intercept = F,standardize = F,lambda = lambda)
  }else if (penalty == "SCAD"){
    current_model <- suppressWarnings(ncvint(as.matrix(logistic_sample[,-1]),logistic_sample[,1],
                                             family="binomial",penalty="SCAD",lambda=lambda,gamma=gamma))
  }else if (penalty == "MCP"){
    current_model <- suppressWarnings(ncvint(as.matrix(logistic_sample[,-1]),logistic_sample[,1],
                                             family="binomial",penalty="MCP",lambda=lambda,gamma=gamma))
  }
  return(t(as.matrix(coef(current_model))))
}
