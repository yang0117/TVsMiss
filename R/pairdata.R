#' Generate paired dataset from original dataset
#'
#' Generate the paired data used for logistics regression
#'
#' @param data1 the original unpaired data
#'
#' @examples
#' n <- 50
#' p <- 8
#' beta <- c(3,0,1.5,0,2,rep(0,p-5))
#' xm <- matrix(rnorm(n*p),ncol = p, nrow = n)
#' y <- xm %*% beta + rnorm(n)
#' colnames(xm) <- paste0("Var_",1:p)
#' pair01 <- pairdata(cbind(y,xm))
#' nrow(pair01) == choose(n,2)
#'
#' @export
pairdata <- function(data1){
  #logistics data generator: get the logistics dataset from original dataset
  #input:original dataset
  #output:logistics dataset

  #covert to logistic regression data
  #define a function to calculate each row of new logistic dataset
  data1 <- as.matrix(data1)
  diff_function <- function(idx_vec,dataset){
    sgn <- dataset[,1][idx_vec[1]] - dataset[,1][idx_vec[2]]
    x_vec <- (dataset[idx_vec[1],-1]-dataset[idx_vec[2],-1])*abs(dataset[idx_vec[1],1]-dataset[idx_vec[2],1])
    res <- c(sgn,x_vec)
    return(res)
  }
  #generate logistic data
  row_idx <- combn(dim(data1)[1],2)
  logistic_sample <- t(apply(row_idx,2,diff_function,dataset=data1))
  logistic_sample <- subset(logistic_sample,logistic_sample[,1] != 0)
  logistic_sample[,1][which(logistic_sample[,1]>0)] <- 1
  logistic_sample[,1][which(logistic_sample[,1]<0)] <- 0
  logistic_sample <- as.matrix(logistic_sample)
  colnames(logistic_sample) <- colnames(data1)

  return(logistic_sample)
}
