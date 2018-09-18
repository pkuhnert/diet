#' Summarises a Node of the Tree 
#' 
#' @description The \code{grab} functions allow you to select nodes of the tree for further 
#' investigation and summary. Once a node is selected, a map of the points is shown
#' with highlighted points of the selection along with a barplot of the prey
#' distribution found at that node. When applied to an object of class \code{bag},
#' bagged distributions of prey are shown. The \code{grabmulti} function produces
#' multiple summaries on successive nodes selected by the user.
#' 
#' @param object tree object either of class \code{dpart} or \code{bag}
#' @param LatID string. Column name of latitude.
#' @param LonID string. Column name of longitude.
#' @param setID optional set identification number.
#' @param node.cols vector of node colours assigned to each prey. 
#' Use output from \code{apc} function. If node colours are not provided, 
#' the trees will not be produced in colour.
#' @param cex numeric. Size of plotting symbols and labels (default: 0.8)
#' @param mapxlim optional map x-limits. If none are specified then the 
#' range of the data are used.
#' @param mapylim optional map y-limits. If none are specified then the 
#' range of the data are used.
#' @param mapcol map land colour. (default: gold3)
#' @param database either 'world' or 'world2' are available for plotting. 
#' Defaults to 'world' if
#' not specified.
#' @param onepage logical. If \code{TRUE}, produces all of the plots on a 
#' single page. (Set to \code{FALSE}
#' if you want publication quality plots appearing on each page).
#' @param display.object tree object to display. This may be different 
#' to \code{object} when bagging is used.
#' @param oob logical. Option when a bagged tree is passed to the 
#' \code{grab.bag} function to 
#' use out of bag estimates.
#' @param n numeric. Maximum number of nodes to investigate in the tree. 
#' Defaults to the number in the tree.
#' @param ylim y-axis limit for barcharts that are produced by the 
#' \code{explore} function
#' @param \dots  arguments to be passed to or from other methods
#' 
#' 
#'                                                        
#' @details The \code{grab.dpart} function can be used on any 
#' tree object of class \code{dpart}.
#'   The \code{grab.bag} function can be used on any tree 
#'   object of class \code{bag}.
#'     The \code{grabmulti.dpart} function can only be used on 
#'     tree objects of class \code{dpart}.
#'       These functions can be invoked explicitly or just by 
#'       calling the \code{grab} 
#'         and \code{grabmulti} functions.
#'         
#' @return   Summary output from the node/s selected in the tree 
#' consisting of the node number,
#' number of observations, number of sets (if the setID was provided), number of
#' predators, number of prey, deviance, expected loss and predicted class.
#' 
#' @references   Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
#' and Stable Isotope Data, Journal of Statistical Software, In Prep. 
#' 
#' @seealso   \code{\link{plot.dpart}}; \code{\link{bagging}}
#' 
#' @examples 
#' # Assigning prey colours for default palette
#'   val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = TRUE)
#'   node.colsY <- val$cols
#'   dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes
#'   
#' # Fitting the classification tree
#'   yft.dp <- dpart(Group ~ Lat + Lon + Year + Quarter + SST  + Length, 
#'                       data = dietPP, weights = W, minsplit = 10,
#'                                         cp = 0.001)
#'                                        
#' # Pruning the tree
#' yft.pr <- prune(yft.dp, se = 1)                   
#' 
#' # Exploring Nodes: This suite of graphics is interactive and therefore has
#' # been commented out. When run, the code will ask you to select a node for
#' # viewing.
#' \dontrun{
#'   # Exploring nodes of the tree - single page
#'   val <- grab(object = yft.pr, LatID = "Lat", LonID = "Lon", setID = "TripSetNo", 
#'               node.cols = node.colsY, cex = 1, mapxlim = c(-125, -75), mapylim = c(0, 30),
#'                           mapcol = "gold3", pos = "topleft")
#'                           
#'  # Exploring nodes of the tree - multiple pages
#'  val <- grab(object = yft.pr, LatID = "Lat", LonID = "Lon", setID = "TripSetNo", 
#'              node.cols = node.colsY, cex = 1, mapxlim = c(-125, -75), mapylim = c(0, 30),
#'                          mapcol = "gold3", pos = "topleft", onepage = FALSE)
#'                          
#'  # Exploring multiple nodes
#'  grabmulti(object = yft.pr, LatID = "Lat", LonID = "Lon", setID = "TripSetNo", 
#'            node.cols = node.colsY, cex = 0.8, mapxlim = c(-125, -75), mapylim = c(0, 30))
#' }
#' @export
grab <- function(object, LatID, LonID, setID = NULL, node.cols = NULL, cex = 0.8, 
                 mapxlim = NULL, mapylim = NULL, mapcol = "gold3", database = 'world', 
                 onepage = TRUE, display.object, oob = FALSE, ylim) UseMethod("grab")


#' @rdname grab
#' @import grid
#' @export
grabmulti.dpart <- function(object, LatID, LonID, setID = NULL, node.cols = NULL, cex = 0.8,
                            mapxlim = NULL, mapylim = NULL, 
                            mapcol = "black", n = nrow(object$frame), database = 'world', ...){
  
  if (!inherits(object, "dpart"))
    stop("Not a dpart object")
  
  if(missing(LatID) | missing(LonID))
    stop("LatID and LonID column names are not specified.\n")
  
  dat <- object$data
  if(is.null(mapxlim))
    mapxlim <- range(dat[,LonID])
  if(is.null(mapylim))
    mapylim <- range(dat[,LatID])
  
  def.par <- par(no.readonly = TRUE)
  
  nodeS <- NULL
  val <- list()
  count <- 0
  while(count < n){
    cat("Click to show the tree\n")
    ll <- grid.locator(unit = "points")
    if (is.null(ll))
      break
    i <- count + 1
    
    val[[i]] <- grab.dpart(object, LatID = LatID, LonID = LonID, setID = setID, node.cols = node.cols, cex = cex,
                           mapxlim = mapxlim, mapylim = mapylim, mapcol = mapcol, database = database, 
                           onepage = TRUE, ...)
    nodeS <- rbind(nodeS, val[[i]]$nodeS)
    count <- count + 1
  }
  
  res <- list(tree = val[[length(val)]]$tree, nodedata = val[[length(val)]]$nodedata,
              nodeS = nodeS)
  
  class(res) <- "grab"
  
  par(def.par)
  
  res
}

#' @rdname grab
#' @export
grabmulti <-function(object, LatID, LonID, setID = NULL, node.cols = NULL, cex = 0.8,
                     mapxlim = NULL, mapylim = NULL, mapcol = "gold3", n = nrow(object$frame), 
                     database = 'world', ...) UseMethod("grabmulti")

#' @rdname grab
#' @export
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

#' @rdname grab
#' @export
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
                                                        aes_string("Lon", "Lat"), pch = 16, size = 2,
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







