grab.dpart <- function(object, LatID, LonID, setID = NULL, node.cols = NULL, cex = 0.8, 
    mapxlim = NULL, mapylim = NULL, mapcol = "black", database = 'world', onepage = TRUE, ...){  
  
  
  if (!inherits(object, "dpart"))
    stop("Not a dpart object")
  
  if(missing(LatID) | missing(LonID))
    stop("LatID and LonID column names are not specified.\n")
  
  
  dat <- object$data
  class(dat) <- c("data", class(dat))
  
  if(is.null(mapxlim))
    mapxlim <- range(dat[,LonID])
  if(is.null(mapylim))
    mapylim <- range(dat[,LatID])
  

  res <- plot.dpart(object, node.cols = node.cols, keep.margins = TRUE)
  
  

val <- rpartco.dpart(object)
node <- as.integer(row.names(object$frame))
temp <- rpart.branch(val$x, val$y, node, branch=1)
nodeLR <- c(temp$nodeL, temp$nodeR)
cat("Click on a node\n")
id <- identify(res$x, res$y, plot = FALSE, n = 1)
nID <- node[id]
cat("Processing information ...\n")


# Extracting observations for plotting on map
subtree <- snip.dpart(object, nID)
nodesub <- as.numeric(row.names(subtree$frame))
val <- rpartco.dpart(subtree)
rn <- row.names(subtree$frame)[subtree$frame$var == "<leaf>"]
pred.where <- names(table(subtree$where))[(1:length(rn))[rn == nID]]
dat.where <- dat[subtree$where == pred.where,]
col.id <- unlist(strsplit(res$labs[id], "\n"))[1]
rect(res$boxes$x1[id], res$boxes$y1[id], res$boxes$x2[id], res$boxes$y2[id],
          col = NA, border = "black", lwd = 2)




# Plot 2: Mapping
m <- mappoints.data(dat[,c(LonID,LatID)], xlim = mapxlim, ylim = mapylim, 
                    database = database) + geom_point(data = dat.where, 
                                                      aes(Lon, Lat), pch = 16, size = 2,
                                                          col = node.cols[col.id]) +
                      ggtitle(paste("Node ", nID, " (n_pred=",
                                    subtree$frame$wt[as.integer(pred.where)],")", sep = ""))
print(m)


# Plot 3: barplot
subpred <- predict(subtree, type = "prob", plot = FALSE)
subpred.where <- subtree$where

# Calculating summary statistics about the node
pred.node <- subpred[subpred.where == paste(pred.where),]
pred.node.m <- apply(pred.node, 2, mean)
node <-  as.integer(row.names(subtree$frame)[as.integer(pred.where)])
nobs <- subtree$frame$n[as.integer(pred.where)]
if(is.null(setID))
  nsets  <- length(unique(subtree$data[subpred.where == as.integer(pred.where),]))
else
  nsets  <- length(unique(subtree$data[subpred.where == as.integer(pred.where),][setID]))
  
npredators <- subtree$frame$wt[as.integer(pred.where)]
nprey <- length(pred.node.m[pred.node.m > 0])
dev <- subtree$frame$dev[as.integer(pred.where)]
loss <- subtree$frame$dev[as.integer(pred.where)]/subtree$frame$wt[as.integer(pred.where)]

pclass <- with(subtree, levels(data$Group)[frame[paste(node),]$yval])
explore(object = subtree, pred = subpred, pred.where = subpred.where, loss = loss,
    node = pred.where, cols = node.cols, showtitle = FALSE, cex = cex)

mtext(side = 3, paste("Diet Composition (D=", round(loss, 3), ")", sep = ""),
    line = 1.5, cex = cex)


res <- list(tree = subtree, nodedata = dat.where,
     nodeS = data.frame(node = node, nobs = nobs, nsets = nsets,
     npredators = npredators, nprey = nprey, dev = dev, loss = loss, pclass = pclass))

class(res) <- "grab"

res

}


