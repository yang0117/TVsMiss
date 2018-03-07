#' plot solution path from the fitted "TVsMiss" object
#'
#' solution path is generated, the x-axis can be either in log or normal scale, the variable names
#' of each predictors can be chosen to show or not
#'
#' @param x fitted "TVsMiss" object
#' @param log If TRUE, x-axis is log scale; if FALSE, x-axis is in normal scale
#' @param label If TRUE, the name of each predictor variable will be showed
#' @param ... graphical parameters to plot
#'
#' @examples
#' n <- 50
#' p <- 8
#' beta <- c(3,0,1.5,0,2,rep(0,p-5))
#' xm <- matrix(rnorm(n*p),ncol = p, nrow = n)
#' y <- xm %*% beta + rnorm(n)
#' colnames(xm) <- paste0("Var_",1:p)
#'
#' fit01 <- tvsmiss(x=xm,y=y)
#' fit01$selection_beta
#' fit01$beta_matrix
#' plot(fit01)
#' plot(fit01,x.log=TRUE,label = FALSE)
#' plot(fit01,x.log=TRUE,label = TRUE)
#'
#' fit04 <- tvsmiss(x=xm,y=y,penalty = "SCAD",method = "BIC")
#' fit04$selection_beta
#' fit04$beta_matrix
#' plot(fit04)
#' plot(fit04,x.log = TRUE)
#' plot(fit04,x.log = TRUE,label = TRUE)
#'
#' @importFrom grDevices hcl
#' @importFrom graphics abline axis mtext par polygon text
#'
#' @export
plot.TVsMiss <- function(x, label=FALSE, log=TRUE, ...){
  if(any(class(x$model) == "glmnet")){
    model.adapt <- x$model
    class(model.adapt) <- "ncvint"
    beta.temp <- as.matrix(model.adapt$beta)
    beta.temp <- rbind(rep(0,ncol(beta.temp)),beta.temp)
    model.adapt$beta <- beta.temp
    model.adapt$penalty.factor <- rep(1,nrow(beta.temp)-1)
    plot_for_ncvreg(model.adapt, log.l=log, label = label,...)
  }else if(any(class(x$model) == "ncvreg")){
    model.adapt <- x$model
    class(model.adapt) <- "ncvint"
    plot_for_ncvreg(model.adapt,log.l = log, label = label,...)
  }else stop("Input obeject is wrong")
}

#ncvreg plot
plot_for_ncvreg <- function(x, alpha=1, log.l=FALSE, shade=FALSE, label=F,...) {
  YY <- if (length(x$penalty.factor)==nrow(x$beta)) coef(x) else coef(x)[-1,,drop=FALSE]
  penalized <- which(x$penalty.factor!=0)
  nonzero <- which(apply(abs(YY), 1, sum)!=0)
  ind <- intersect(penalized, nonzero)
  Y <- YY[ind, , drop=FALSE]
  p <- nrow(Y)
  l <- x$lambda

  #print(Y)

  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)
  plot.args <- list(x=l, y=1:length(l), ylim=range(Y), xlab=xlab, ylab="", type="n", xlim=c(rev(range(c(l,min(l)-abs(max(l)-min(l))*0.1)))), las=1)
  new.args <- list(...)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("matplot", plot.args)
  if (!is.element("ylab", names(new.args))) mtext(expression(hat(beta)), side=2, cex=par("cex"), line=3, las=1)
  if (shade & !is.null(x$convex.min)) {
    l1 <- l[x$convex.min]
    l2 <- min(l)
    polygon(x=c(l1,l2,l2,l1),y=c(plot.args$ylim[1],plot.args$ylim[1],plot.args$ylim[2],plot.args$ylim[2]),col="gray85",border=FALSE)
  }
  #axis(3, col = "gold", lty = 2, lwd = 0.5)

  index = l
  atdf = pretty(index)
  df = as.vector(apply(x$beta,2, FUN=function(s) sum(s[-1]!= 0)))
  approx.f = 0
  prettydf = approx(x = index, y = df, xout = atdf, rule = 2,
                    method = "constant", f = approx.f)$y

  axis(3, at = atdf, labels = prettydf, tcl = NA)

  cols <- hcl(h=seq(15, 375, len=max(4, p+1)), l=60, c=150, alpha=alpha)
  cols <- if (p==2) cols[c(1,3)] else cols[1:p]
  line.args <- list(col=cols, lwd=1+2*exp(-p/20), lty=1)
  if (length(new.args)) line.args[names(new.args)] <- new.args
  line.args$x <- l
  line.args$y <- t(Y)
  do.call("matlines",line.args)

  if (label) {
    beta <- x$beta

    # if(log.l){
    # xpos <- min(log(index))
    # }else{
    # xpos <- min(index)
    # }

    xpos <- min(index)
    xpos <- rep(xpos, (nrow(beta)-1))
    ypos = beta[2:nrow(beta), ncol(beta)]

    pos = 4

    text(xpos, ypos, rownames(beta)[2:length(rownames(beta))], cex = 0.5, pos = pos)
  }
  abline(h=0)
}
