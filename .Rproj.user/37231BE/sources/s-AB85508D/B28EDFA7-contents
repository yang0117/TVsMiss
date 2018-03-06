# Function 01

# three steps:
# 1. remove missing and pairwise each observations
# 2. use penalty to get lambda path and corresponding beta matrix

# mispairwise <- function(x,y){
#
#   this.call = match.call()
#
#   if (class(x) != "matrix") {
#     tmp <- try(x <- model.matrix(~0+., data=x), silent=TRUE)
#     if (class(tmp)[1] == "try-error") stop("x must be a matrix or able to be coerced to a matrix")
#   }
#   if (storage.mode(x)=="integer") storage.mode(x) <- "double"
#   if (class(y) != "numeric") {
#     tmp <- try(y <- as.numeric(y), silent=TRUE)
#     if (class(tmp)[1] == "try-error") stop("y must numeric or able to be coerced to numeric")
#   }
#
#   if (length(y) != nrow(x)) stop("x and y do not have the same number of observations")
#
#   #step 1: remove missing
#   colnames(x) <- if (is.null(colnames(x))) paste("V",1:ncol(x),sep="") else colnames(x)
#
#   data <- cbind(y,x)
#   complete_idx <- c(1:length(y))[complete.cases(data)]
#   y_complete <- y[complete_idx]
#   x_complete <- x[complete_idx,]
#
#   #generate logistic data
#   sample_no_missing <- cbind(y_complete,x_complete)
#   logistic_sample <- pairdata(sample_no_missing)
#   current_model <- glm(logistic_sample[,1]~logistic_sample[,-1]-1, family = "binomial")
#   selection_beta <- coef(current_model)
#
#   deviance.ratio = 1 + loglikelihood(c(0,selection_beta),sample_no_missing)/log(2)
#
#   res <- list(ls=logistic_sample,c_idx=complete_idx,model=current_model,
#               selection_beta =  selection_beta,
#               null.deviance = 2*log(2),
#               deviance.ratio=deviance.ratio,
#               call = this.call)
#   return(res)
# }

