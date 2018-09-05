link.bag <- function(x, object, LatID, LonID, mapxlim = NULL, mapylim = NULL,
   plot = TRUE,  oob = FALSE, mfrow = c(2,2), orderN = FALSE){

  # links the predictions to the terminal nodes of a tree model
  options(warn = -1)
  if (!inherits(x, "bag")) 
        stop("Not a bagged object")
  
  dat <- object$data
  if(plot){ 
    if(is.null(mapxlim))
      mapxlim <- range(dat[,LonID])
    if(is.null(mapylim))
      mapylim <- range(dat[,LatID])
  } 

  nBaggs <- length(x$baggs)
  nodenms <- row.names(object$frame)[object$frame$var == "<leaf>"]
 

  bpred <- list()
  for(i in 1:nBaggs){
      bpred[[i]] <- data.frame(matrix(0, nrow = length(names(table(object$where))),
                ncol = ncol(x$pred[[i]])))
#         bpred[[i]] <- data.frame(matrix(0, nrow = length(nodenms),
#                    ncol = ncol(x$pred[[i]])))
    row.names(bpred[[i]]) <- names(table(object$where))
      #   row.names(bpred[[i]]) <- nodenms
      names(bpred[[i]]) <- names(x$pred[[i]])
  }
  options(warn = -1)
  for(i in 1:nBaggs){
      if(oob){  # need to ensure all nodes are represented in the list 
       where <- rpart:::pred.rpart(object, rpart.matrix(x$data[paste(x$oob[[i]]),]))
         tmp <- apply(x$pred.oob[[i]], 2, function(x, wh) tapply(x, wh, mean), where)
         n.nms <- names(table(where))
         bpred[[i]][n.nms,] <- tmp
      }
      else{
         where <- object$where
         bpred[[i]] <- apply(x$pred[[i]], 2, function(x, wh) tapply(x, wh, mean), where)
      }
  }

  bpred.sum <- do.call(abind, c(bpred, along = 3))
  bpred.m <- apply(bpred.sum, c(1,2), mean)
  bpred.v <- apply(bpred.sum, c(1,2), var)
  bpred.lci <- apply(bpred.sum, c(1,2), function(x) quantile(x, 0.025))
  bpred.uci <- apply(bpred.sum, c(1,2), function(x) quantile(x, 0.975))
  # ordering the node outputs
  if(orderN){
     bpred.m <- bpred.m[order(as.numeric(nodenms)),]
     bpred.v <- bpred.v[order(as.numeric(nodenms)),]
     bpred.lci <- bpred.lci[order(as.numeric(nodenms)),]
     bpred.uci <- bpred.uci[order(as.numeric(nodenms)),]
     nodenms <- as.character(sort(as.numeric(nodenms)))
  }

  if(plot){
  if(dev.cur() == 1) windows(8,8, record = TRUE)
  par(mfrow=mfrow)
  nc <- ncol(bpred.m)

  ylimit <- c(-0.05,1.05)

  options(warn = -1)
  for(i in 1:nrow(bpred.m)){
      plot(1:nc, bpred.m[i,], pch = 15, ylim = ylimit, axes = FALSE,
         xlab = "", ylab = "Proportion")
      box()
      axis(side = 2)
      mtext(side = 1, text = names(data.frame(bpred.m)), at = 1:nc, line = 1, 
        adj = 1, cex = 0.8, las = 2)
      arrows(1:nc, bpred.lci[i,], 1:nc, bpred.uci[i,], length = 0.1, angle = 90,
        code = 3)  
      title(main = paste("Node ", nodenms[i], sep = ""))
  }
  par(mfrow=c(1,1))
  }
 
  row.names(bpred.m) <- nodenms
  row.names(bpred.v) <- nodenms
  row.names(bpred.lci) <- nodenms
  row.names(bpred.uci) <- nodenms
  options(warn = 0)
  
  val <- list(m = bpred.m, v = bpred.v, lci95 = bpred.lci, uci95 = bpred.uci)
  
  class(val) <- c("link")
  
  val

} 

