#' RPART functions
#' 
#' @description These have been exported from the rpart package to
#' avoid errors in the check
#' 
#' 

#' @rdname rpart-Internal
rpart.class <- function (y, offset, parms, wt) 
{
  if (!is.null(offset)) 
    stop("No offset allowed in classification models")
  fy <- as.factor(y)
  y <- as.integer(fy)
  numclass <- max(y[!is.na(y)])
  counts <- tapply(wt, factor(y, levels = 1:numclass), sum)
  counts <- ifelse(is.na(counts), 0, counts)
  if (missing(parms) || is.null(parms)) 
    parms <- list(prior = counts/sum(counts), loss = matrix(rep(1, 
                                                                numclass^2) - diag(numclass), numclass), split = 1)
  else if (is.list(parms)) {
    if (is.null(names(parms))) 
      stop("The parms list must have names")
    temp <- pmatch(names(parms), c("prior", "loss", "split"), 
                   0L)
    if (any(temp == 0L)) 
      stop(gettextf("'parms' component not matched: %s", 
                    names(parms)[temp == 0L]), domain = NA)
    names(parms) <- c("prior", "loss", "split")[temp]
    if (is.null(parms$prior)) 
      temp <- c(counts/sum(counts))
    else {
      temp <- parms$prior
      if (sum(temp) != 1) 
        stop("Priors must sum to 1")
      if (any(temp < 0)) 
        stop("Priors must be >= 0")
      if (length(temp) != numclass) 
        stop("Wrong length for priors")
    }
    if (is.null(parms$loss)) 
      temp2 <- 1 - diag(numclass)
    else {
      temp2 <- parms$loss
      if (length(temp2) != numclass^2) 
        stop("Wrong length for loss matrix")
      temp2 <- matrix(temp2, ncol = numclass)
      if (any(diag(temp2) != 0)) 
        stop("Loss matrix must have zero on diagonals")
      if (any(temp2 < 0)) 
        stop("Loss matrix cannot have negative elements")
      if (any(rowSums(temp2) == 0)) 
        stop("Loss matrix has a row of zeros")
    }
    if (is.null(parms$split)) 
      temp3 <- 1L
    else {
      temp3 <- pmatch(parms$split, c("gini", "information"))
      if (is.null(temp3)) 
        stop("Invalid splitting rule")
    }
    parms <- list(prior = temp, loss = matrix(temp2, numclass), 
                  split = temp3)
  }
  else stop("Parameter argument must be a list")
  list(y = y, parms = parms, numresp = numclass + 2L, counts = counts, 
       ylevels = levels(fy), numy = 1L, print = function(yval, 
                                                         ylevel, digits) {
         temp <- if (is.null(ylevel)) as.character(yval[, 
                                                        1L]) else ylevel[yval[, 1L]]
         nclass <- (ncol(yval) - 2L)/2L
         yprob <- if (nclass < 5L) format(yval[, 1L + nclass + 
                                                 1L:nclass], digits = digits, nsmall = digits) else formatg(yval[, 
                                                                                                                 1L + nclass + 1L:nclass], digits = 2L)
         if (!is.matrix(yprob)) yprob <- matrix(yprob, nrow = 1L)
         temp <- paste0(temp, " (", yprob[, 1L])
         for (i in 2L:ncol(yprob)) temp <- paste(temp, yprob[, 
                                                             i], sep = " ")
         temp <- paste0(temp, ")")
         temp
       }, summary = function(yval, dev, wt, ylevel, digits) {
         nclass <- (ncol(yval) - 2L)/2L
         group <- yval[, 1L]
         counts <- yval[, 1L + (1L:nclass)]
         yprob <- yval[, 1L + nclass + 1L:nclass]
         nodeprob <- yval[, 2L * nclass + 2L]
         if (!is.null(ylevel)) group <- ylevel[group]
         temp1 <- formatg(counts, format = "%5g")
         temp2 <- formatg(yprob, format = "%5.3f")
         if (nclass > 1) {
           temp1 <- apply(matrix(temp1, ncol = nclass), 
                          1L, paste, collapse = " ")
           temp2 <- apply(matrix(temp2, ncol = nclass), 
                          1L, paste, collapse = " ")
         }
         dev <- dev/(wt[1L] * nodeprob)
         paste0("  predicted class=", format(group, justify = "left"), 
                "  expected loss=", formatg(dev, digits), "  P(node) =", 
                formatg(nodeprob, digits), "\n", "    class counts: ", 
                temp1, "\n", "   probabilities: ", temp2)
       }, text = function(yval, dev, wt, ylevel, digits, n, 
                          use.n) {
         nclass <- (ncol(yval) - 2L)/2L
         group <- yval[, 1L]
         counts <- yval[, 1L + (1L:nclass)]
         if (!is.null(ylevel)) group <- ylevel[group]
         temp1 <- formatg(counts, digits)
         if (nclass > 1L) temp1 <- apply(matrix(temp1, ncol = nclass), 
                                         1L, paste, collapse = "/")
         if (use.n) paste0(format(group, justify = "left"), 
                           "\n", temp1) else format(group, justify = "left")
       })
}

#' @rdname rpart-Internal
rpart.anova <- function (y, offset, parms, wt) 
{
  if (!is.null(offset)) 
    y <- y - offset
  list(y = y, parms = NULL, numresp = 1L, numy = 1L, summary = function(yval, 
                                                                        dev, wt, ylevel, digits) {
    paste0("  mean=", formatg(yval, digits), ", MSE=", formatg(dev/wt, 
                                                               digits))
  }, text = function(yval, dev, wt, ylevel, digits, n, use.n) {
    if (use.n) paste0(formatg(yval, digits), "\nn=", n) else formatg(yval, 
                                                                     digits)
  })
}

#' @rdname rpart-Internal
rpart.poisson <- function (y, offset, parms, wt) 
{
  if (is.matrix(y)) {
    if (ncol(y) != 2L) 
      stop("response must be a 2 column matrix or a vector")
    if (!is.null(offset)) 
      y[, 1L] <- y[, 1L] * exp(offset)
  }
  else {
    if (is.null(offset)) 
      y <- cbind(1, y)
    else y <- cbind(exp(offset), y)
  }
  if (any(y[, 1L] <= 0)) 
    stop("Observation time must be > 0")
  if (any(y[, 2L] < 0)) 
    stop("Number of events must be >= 0")
  if (missing(parms)) 
    parms <- c(shrink = 1L, method = 1L)
  else {
    parms <- as.list(parms)
    if (is.null(names(parms))) 
      stop("You must input a named list for parms")
    parmsNames <- c("method", "shrink")
    indx <- pmatch(names(parms), parmsNames, 0L)
    if (any(indx == 0L)) 
      stop(gettextf("'parms' component not matched: %s", 
                    names(parms)[indx == 0L]), domain = NA)
    else names(parms) <- parmsNames[indx]
    if (is.null(parms$method)) 
      method <- 1L
    else method <- pmatch(parms$method, c("deviance", "sqrt"))
    if (is.null(method)) 
      stop("Invalid error method for Poisson")
    if (is.null(parms$shrink)) 
      shrink <- 2L - method
    else shrink <- parms$shrink
    if (!is.numeric(shrink) || shrink < 0L) 
      stop("Invalid shrinkage value")
    parms <- c(shrink = shrink, method = method)
  }
  list(y = y, parms = parms, numresp = 2L, numy = 2L, summary = function(yval, 
                                                                         dev, wt, ylevel, digits) {
    paste0("  events=", formatg(yval[, 2L]), ",  estimated rate=", 
           formatg(yval[, 1L], digits), " , mean deviance=", 
           formatg(dev/wt, digits))
  }, text = function(yval, dev, wt, ylevel, digits, n, use.n) {
    if (!is.matrix(yval)) yval <- matrix(yval, nrow = 1L)
    if (use.n) paste0(formatg(yval[, 1L], digits), "\n", 
                      formatg(yval[, 2L]), "/", n) else paste(formatg(yval[, 
                                                                           1L], digits))
  })
}

#' @rdname rpart-Internal
rpartco <- function (tree, parms) 
{
  if (missing(parms)) {
    pn <- paste0("device", dev.cur())
    if (!exists(pn, envir = rpart:::rpart_env, inherits = FALSE)) 
      stop("no information available on parameters from previous call to plot()")
    parms <- get(pn, envir = rpart:::rpart_env, inherits = FALSE)
  }
  frame <- tree$frame
  node <- as.numeric(row.names(frame))
  depth <- tree.depth(node)
  is.leaf <- (frame$var == "<leaf>")
  if (length(parms)) {
    uniform <- parms$uniform
    nspace <- parms$nspace
    minbranch <- parms$minbranch
  }
  else {
    uniform <- FALSE
    nspace <- -1
    minbranch <- 0.3
  }
  if (uniform) 
    y <- (1 + max(depth) - depth)/max(depth, 4L)
  else {
    y <- dev <- frame$dev
    temp <- split(seq(node), depth)
    parent <- match(node%/%2L, node)
    sibling <- match(ifelse(node%%2L, node - 1L, node + 1L), 
                     node)
    for (i in temp[-1L]) {
      temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
      y[i] <- y[parent[i]] - temp2
    }
    fudge <- minbranch * diff(range(y))/max(depth)
    for (i in temp[-1L]) {
      temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
      haskids <- !(is.leaf[i] & is.leaf[sibling[i]])
      y[i] <- y[parent[i]] - ifelse(temp2 <= fudge & haskids, 
                                    fudge, temp2)
    }
    y <- y/(max(y))
  }
  x <- double(length(node))
  x[is.leaf] <- seq(sum(is.leaf))
  left.child <- match(node * 2L, node)
  right.child <- match(node * 2L + 1L, node)
  temp <- split(seq(node)[!is.leaf], depth[!is.leaf])
  for (i in rev(temp)) x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])
  if (nspace < 0) 
    return(list(x = x, y = y))
  compress <- function(x, me, depth) {
    lson <- me + 1L
    if (is.leaf[lson]) 
      left <- list(left = x[lson], right = x[lson], depth = depth + 
                     1L, sons = lson)
    else {
      left <- compress(x, me + 1L, depth + 1L)
      x <- left$x
    }
    rson <- me + 1L + length(left$sons)
    if (is.leaf[rson]) 
      right <- list(left = x[rson], right = x[rson], depth = depth + 
                      1L, sons = rson)
    else {
      right <- compress(x, rson, depth + 1L)
      x <- right$x
    }
    maxd <- max(left$depth, right$depth) - depth
    mind <- min(left$depth, right$depth) - depth
    slide <- min(right$left[1L:mind] - left$right[1L:mind]) - 
      1L
    if (slide > 0) {
      x[right$sons] <- x[right$sons] - slide
      x[me] <- (x[right$sons[1L]] + x[left$sons[1L]])/2
    }
    else slide <- 0
    if (left$depth > right$depth) {
      templ <- left$left
      tempr <- left$right
      tempr[1L:mind] <- pmax(tempr[1L:mind], right$right - 
                               slide)
    }
    else {
      templ <- right$left - slide
      tempr <- right$right - slide
      templ[1L:mind] <- pmin(templ[1L:mind], left$left)
    }
    list(x = x, left = c(x[me] - nspace * (x[me] - x[lson]), 
                         templ), right = c(x[me] - nspace * (x[me] - x[rson]), 
                                           tempr), depth = maxd + depth, sons = c(me, left$sons, 
                                                                                  right$sons))
  }
  x <- compress(x, 1L, 1L)$x
  list(x = x, y = y)
}

formatg <- function (x, digits = getOption("digits"), format = paste0("%.", 
                                                           digits, "g")) 
{
  if (!is.numeric(x)) 
    stop("'x' must be a numeric vector")
  temp <- sprintf(format, x)
  if (is.matrix(x)) 
    matrix(temp, nrow = nrow(x))
  else temp
}

#' @rdname rpart-Internal
plot.rpart <- function (x, uniform = FALSE, branch = 1, compress = FALSE, nspace, 
                        margin = 0, minbranch = 0.3, ...) 
{
  if (!inherits(x, "rpart")) 
    stop("Not a legitimate \"rpart\" object")
  if (nrow(x$frame) <= 1L) 
    stop("fit is not a tree, just a root")
  if (compress & missing(nspace)) 
    nspace <- branch
  if (!compress) 
    nspace <- -1L
  parms <- list(uniform = uniform, branch = branch, nspace = nspace, 
                minbranch = minbranch)
  temp <- rpartco(x, parms)
  xx <- temp$x
  yy <- temp$y
  temp1 <- range(xx) + diff(range(xx)) * c(-margin, margin)
  temp2 <- range(yy) + diff(range(yy)) * c(-margin, margin)
  plot(temp1, temp2, type = "n", axes = FALSE, xlab = "", ylab = "", 
       ...)
  assign(paste0("device", dev.cur()), parms, envir = rpart:::rpart_env)
  node <- as.numeric(row.names(x$frame))
  temp <- rpart.branch(xx, yy, node, branch)
  if (branch > 0) 
    text(xx[1L], yy[1L], "|")
  lines(c(temp$x), c(temp$y))
  invisible(list(x = xx, y = yy))
}


