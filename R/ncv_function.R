lamNames <- function(l) {
  if (length(l) > 1) {
    d <- ceiling(-log10(-max(diff(l))))
    d <- min(max(d,4), 10)
  } else {
    d <- 4
  }
  formatC(l, format="f", digits=d)
}

convexMin <- function(b, X, penalty, gamma, l2, family, penalty.factor, a, Delta=NULL) {
  n <- nrow(X)
  p <- ncol(X)
  l <- ncol(b)

  if (penalty=="MCP") {
    k <- 1/gamma
  } else if (penalty=="SCAD") {
    k <- 1/(gamma-1)
  } else if (penalty=="lasso") {
    return(NULL)
  }
  if (l==0) return(NULL)

  val <- NULL
  for (i in 1:l) {
    A1 <- if (i==1) rep(1,p) else b[,i]==0
    if (i==l) {
      L2 <- l2[i]
      U <- A1
    } else {
      A2 <- b[,i+1]==0
      U <- A1&A2
      L2 <- l2[i+1]
    }
    if (sum(!U)==0) next
    Xu <- X[,!U]
    p.. <- k*(penalty.factor[!U]!=0) - L2*penalty.factor[!U]
    if (family=="gaussian") {
      if (any(A1!=A2)) {
        eigen.min <- min(eigen(crossprod(Xu)/n - diag(p.., length(p..), length(p..)))$values)
      }
    } else if (family=="binomial") {
      if (i==l) eta <- a[i] + X%*%b[,i]
      else eta <- a[i+1] + X%*%b[,i+1]
      pi. <- exp(eta)/(1+exp(eta))
      w <- as.numeric(pi.*(1-pi.))
      w[eta > log(.9999/.0001)] <- .0001
      w[eta < log(.0001/.9999)] <- .0001
      Xu <- sqrt(w) * cbind(1,Xu)
      xwxn <- crossprod(Xu)/n
      eigen.min <- min(eigen(xwxn-diag(c(0,diag(xwxn)[-1]*p..)))$values)
    } else if (family=="poisson") {
      if (i==l) eta <- a[i] + X%*%b[,i]
      else eta <- a[i+1] + X%*%b[,i+1]
      mu <- exp(eta)
      w <- as.numeric(mu)
      Xu <- sqrt(w) * cbind(1,Xu)
      xwxn <- crossprod(Xu)/n
      eigen.min <- min(eigen(xwxn-diag(c(0,diag(xwxn)[-1]*p..)))$values)
    } else if (family=="cox") {
      eta <- if (i==l) X%*%b[,i] else X%*%b[,i+1]
      haz <- drop(exp(eta))
      rsk <- rev(cumsum(rev(haz)))
      h <- haz*cumsum(Delta/rsk)
      xwxn <- crossprod(sqrt(h) * Xu)/n
      eigen.min <- min(eigen(xwxn-diag(diag(xwxn)*p.., nrow(xwxn), ncol(xwxn)))$values)
    }

    if (eigen.min < 0) {
      val <- i
      break
    }
  }
  val
}

coef.ncvint <- function(object, lambda, which=1:length(object$lambda), drop=TRUE, ...) {
  if (!missing(lambda)) {
    ind <- approx(object$lambda,seq(object$lambda),lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    beta <- (1-w)*object$beta[,l,drop=FALSE] + w*object$beta[,r,drop=FALSE]
    colnames(beta) <- lamNames(lambda)
  }
  else beta <- object$beta[, which, drop=FALSE]
  if (drop) return(drop(beta)) else return(beta)
}

setupLambda <- function(X, y, family, alpha, lambda.min, nlambda, penalty.factor) {
  n <- nrow(X)
  p <- ncol(X)

  ## Determine lambda.max
  ind <- which(penalty.factor!=0)
  if (length(ind)!=p) {
    fit <- glm(y~X[, -ind], family=family)
  } else {
    fit <- glm(y~1, family=family)
  }
  if (family=="gaussian") {
    zmax <- .Call("maxprod", X, fit$residuals, ind, penalty.factor) / n
  } else {
    zmax <- .Call("maxprod", X, residuals(fit, "working") * fit$weights, ind, penalty.factor) / n
  }
  lambda.max <- zmax/alpha

  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
  } else {
    lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
  }

  if (length(ind)!=p) lambda[1] <- lambda[1] * 1.000001
  lambda
}
