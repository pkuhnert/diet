textG.dpart <-
  function (x, splits = TRUE, which = 4, label = "yval", FUN = text,
            all.leaves = FALSE, pretty = NULL, digits = getOption("digits") - 
              2, tadj = 0.65,  use.n = FALSE, bars = TRUE,
            xadj = 1, yadj = 1, bord = FALSE, node.cols = NULL, pos = NULL,
            ...) 
  {
  
    if (!inherits(x, "rpart")) 
      stop("Not legitimate rpart")
  #  if (!is.null(x$frame$splits)) 
  #    x <- rpconvert(x)  # commented this out as compiler didn't like this
    frame <- x$frame
    col <- names(frame)
    method <- x$method
    ylevels <- attr(x, "ylevels")
    if (!is.null(ylevels <- attr(x, "ylevels"))) 
      col <- c(col, ylevels)
    if (is.na(match(label, col))) 
      stop("Label must be a column label of the frame component of the tree")
    cxy <- par("cxy")
    if (!is.null(srt <- list(...)$srt) && srt == 90) 
      cxy <- rev(cxy)
    
    xy <- rpartco.dpart(x)
    node <- as.numeric(row.names(x$frame))
    is.left <- (node%%2 == 0)
    node.left <- node[is.left]
    parent <- match(node.left/2, node)
    
    if (splits) {
      left.child <- match(2 * node, node)
      right.child <- match(node * 2 + 1, node)
      rows <- labels(x, pretty = pretty)
      if (which == 1) 
        FUN(xy$x, xy$y + tadj * cxy[2], rows[left.child], 
            ...)
      else {
        if (which == 2 | which == 4) 
          FUN(xy$x, xy$y + tadj * cxy[2], rows[left.child], 
              pos = 2, ...)
        if (which == 3 | which == 4) 
          FUN(xy$x, xy$y + tadj * cxy[2], rows[right.child], 
              pos = 4, ...)
      }
    }
    leaves <- if (all.leaves) 
      rep(TRUE, nrow(frame))
    else frame$var == "<leaf>"
    
    if (is.null(frame$yval2))
      stat <- x$functions$text(yval = frame$yval[leaves],
                               dev = frame$dev[leaves], wt = frame$wt[leaves],
                               ylevel = ylevels, digits = digits, n = frame$n[leaves],
                               use.n = use.n)
    else stat <- x$functions$text(yval = frame$yval2[leaves,
                                                     ], dev = frame$dev[leaves], wt = frame$wt[leaves],
                                  ylevel = ylevels, digits = digits, n = frame$n[leaves],
                                  use.n = use.n)
    
    if(is.null(node.cols)){
      stat.cols <- as.factor(stat)
      levels(stat.cols) <- rainbow(length(levels(stat.cols)))
    }
    else{
      stat.cols <- as.factor(unlist(lapply(strsplit(stat, " "), function(x) x[1])))
      levels(stat.cols) <-  node.cols[match(levels(stat.cols), names(node.cols))]
    }
    stat <- as.factor(stat)
    
    # Plotting
    text.adj <- yadj * diff(range(xy$y))/50
    FUN(xy$x[leaves], xy$y[leaves] - tadj * cxy[2] - text.adj,
        stat, adj = 0.5, col = "black", ...)
    
    points(xy$x[leaves], xy$y[leaves], pch = 16, cex = 3 *
      par()$cex, col = as.vector(stat.cols))
    
    if(is.null(pos)){
      rx <- range(xy$x)
      ry <- range(xy$y)
      x.leg <- min(xy$x) - 0 * rx
      y.leg <- max(xy$y) + 0 * ry
    }
    
       if(is.null(node.cols)){
           if(is.null(pos))
             legend(x.leg, y.leg, levels(stat), col = levels(stat.cols),
                pch = 16, bty = "n", ...)
           else
             legend(pos, levels(stat), col = levels(stat.cols),
                pch = 16, bty = "n", ...)
        }
        else{
           if(is.null(pos))
              legend(x.leg, y.leg, names(node.cols), col = node.cols,
                 pch = 16, bty = "n", ...)
           else
              legend(pos, names(node.cols), col = node.cols,
                 pch = 16, bty = "n", ...)
        }
    
    invisible()
  }
