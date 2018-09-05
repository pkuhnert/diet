grab.bag <- function(object, LatID, LonID, setID =  NULL, node.cols = NULL, cex = 0.8,
                       mapxlim = NULL, mapylim = NULL, mapcol = "gold3", database = 'world',  
                     onepage = TRUE, display.object, oob = FALSE, ylim){
    

if (!inherits(object, "bag"))
    stop("Not a bagged object. \n")

if(missing(display.object))
  stop("display.object not set. \n")

if(missing(LatID) | missing(LonID))
  stop("LatID and LonID column names are not specified.\n")

def.par <- par(no.readonly = TRUE)
dat <- object$data

if(is.null(mapxlim))
  mapxlim <- range(dat[,LonID])
if(is.null(mapylim))
  mapylim <- range(dat[,LatID])



#if(onepage){
#  layout(rbind( c(1, 1, 0, 2, 2),
#                c(1, 1, 0, 2, 2), 
#                c(1, 1, 0, 3, 3),
#                c(1, 1, 0, 3, 3)), widths = c(1.5, 1.5, lcm(1.5), 1.5, 1), 
#                heights = c(2, 1, 2, 1), 
#         respect = TRUE)
  
#  par(cex = 0.75)
  # Plot 1: Tree
#  plotG.dpart(display.object, node.cols = node.cols, cex = cex, pos = "topleft")
  
#}
#else{
  # Plot 1: Tree
#  plot.dpart(display.object, node.cols = node.cols, keep.margins = TRUE)
#  
#}

res <- plot.dpart(display.object, node.cols = node.cols, keep.margins = TRUE)

val <- rpartco.dpart(display.object)
node <- as.integer(row.names(display.object$frame))
temp <- rpart.branch(val$x, val$y, node, branch=1)
nodeLR <- c(temp$nodeL, temp$nodeR)
cat("Click on a node\n")
id <- identify(res$x, res$y, plot = FALSE, n = 1)
nID <- node[id]
cat("Processing information ...\n")

# Extracting observations for plotting on map
subtree <- snip.dpart(display.object, nID)
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



loss <- subtree$frame$dev[as.integer(pred.where)]/subtree$frame$wt[as.integer(pred.where)]

# barplot
subpred <- predict(subtree, type = "prob", plot = FALSE)
subpred.where <- subtree$where
explore(object = subtree, pred = subpred, pred.where = subpred.where, loss,
    node = pred.where, cols = node.cols, showtitle = FALSE, labels = FALSE, cex = cex, ylim = ylim)


# bootstrap estimates
if(missing(ylim))
  ylim <- c(-0.05,1.05)
bag.map <- link(object, subtree, LatID, LonID,  mapxlim, mapylim, oob = oob, plot = FALSE)
options(warn = -1)

nodestats <- explore.bag(bag.map, nID, cols = node.cols,
     showtitle = FALSE, axis.side = 4, cex = cex, ylim = ylim)
legend("topleft", legend = "Bootstrapped proportions", bty = "n")




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
pclass <- with(subtree, levels(data$Group)[frame[paste(node),]$yval])

par(def.par)

res <- list(tree = subtree, nodedata = dat.where, bag.map = bag.map,
     nodeS = data.frame(node = node, nobs = nobs, nsets = nsets,
     npredators = npredators, nprey = nprey, dev = dev, loss = loss, pclass = pclass),
     nodestats = nodestats)

class(res) <- "grab"

options(warn = 0)

res

}


