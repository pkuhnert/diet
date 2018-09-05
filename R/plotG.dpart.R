plotG.dpart <- function(x, node.cols = NULL, pos = NULL, ...){
  

#  plot(x, keep.margins = TRUE, ...)
  rpart:::plot.rpart(x, ...)
  textG.dpart(x, xpd = NA, pretty = TRUE, splits = TRUE, node.cols = node.cols, 
              pos = pos,  ...)
  
  
  val <- rpartco.dpart(x)
  ff <- x$frame
  nodes <- as.numeric(row.names(ff))
  temp <- rpart.branch(val$x, val$y, nodes, branch=1)
  x1 <- c(temp$x)[seq(1, length(c(temp$x)), by = 5)]
  y1 <- c(temp$y)[seq(1, length(c(temp$y)), by = 5)]
  x2 <- c(temp$x)[seq(4, length(c(temp$x)), by = 5)]
  y2 <- c(temp$y)[seq(4, length(c(temp$y)), by = 5)]
  xvec <- c(x1, x2)
  yvec <- c(y1, y2)
  nodeLR <- c(temp$nodeL, temp$nodeR)
  
  id <- seq_along(xvec)
  nID <- nodeLR[id]
  rn <- row.names(ff)[-1]
  ff.var <- as.vector(ff$var[-1])[match(nID, rn)]
  
  text(xvec[id][ff.var == "<leaf>"], yvec[id][ff.var == "<leaf>"],
       nID[ff.var == "<leaf>"], cex = 0.7, col = "black")
  
 
  
}
