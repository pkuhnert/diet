rpconvert <- function (x) 
{
  if (!inherits(x, "rpart")) 
    stop("x does not appear to be an rpart object")
  ff <- x$frame
  if (is.null(ff$splits)) {
    warning("x not converted")
    return(x)
  }
  ff$splits <- NULL
  ff$wt <- ff$n
  xlev <- attr(x, "xlevels")
  if (length(xlev) > 0) {
    zz <- as.numeric(names(xlev))
    names(xlev) <- attr(x$terms, "term.labels")[zz]
    attr(x, "xlevels") <- xlev
  }
  if (x$method == "class") {
    temp <- cbind(ff$yval, ff$yval2, ff$yprob)
    dimnames(temp) <- NULL
    ff$yval2 <- temp
    ff$yprob <- NULL
    x$frame <- ff
    temp <- rpart:::rpart.class(c(1, 1, 2, 2), NULL, wt = c(1, 1, 
                                                    1, 1))
    x$functions <- list(summary = temp$summary, print = temp$print, 
                        text = temp$text)
  }
  else if (x$method == "anova") {
    x$frame <- ff
    temp <- rpart:::rpart.anova(1:5, NULL, wt = rep(1, 5))
    x$functions <- list(summary = temp$summary, text = temp$text)
  }
  else {
    ff$yval2 <- cbind(ff$yval, ff$yval2)
    x$frame <- ff
    temp <- rpart:::rpart.poisson(1:5, NULL, wt = rep(1, 5))
    x$functions <- list(summary = temp$summary, text = temp$text)
  }
  class(x) <- "rpart"
  x
}