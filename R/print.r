#' Print from a dpart object
#' 
#' @description Either prints a \code{dpart} object or prints output from a node
#' selected from a \code{rpart} tree.
#' 
#' @usage 
#' \method{print}{dpart}(x, setID = NULL, digits=getOption("digits"), 
#'                       file = "diet_tree_summary.csv", ...)
#' \method{print}{grab}(x, ...)
#' 
#' @param x fitted model object of class \code{dpart} OR object of class \code{grab}. 
#' @param setID  the set identification number used to summarise the output. Defaults to NULL.
#' @param digits  the number of digits of numbers to print
#' @param file file name where to output summary results. (Output too big to print to screen).
#' @param \dots  arguments to be passed to or from other methods.
#' 
#' @param details This function is a method for the generic function \code{print} for class "dpart". 
#' It can be invoked by calling print for an object of the appropriate class, or
#' directly by calling print.dpart. The loss is computed as the deviance at a node divided by the 
#' number of predators appearing at a node and ranges between 0 and 1. The loss is a measure of the 
#' diversity of prey eaten by predators at the node, where values near 0 indicate
#' low diversity and values near 1 indicate high diversity.
#' 
#' @details \code{print.dpart} Summary of the tree is provided as a csv file with the following column headings:
#' \item{node}{node number}
#' \item{nobs}{number of observations in each node}
#' \item{nsets}{number of sets in each node if a setID was provided otherwise, NA}
#' \item{npredator}{number of predators in each node}
#' \item{nprey}{number of prey in each node}
#' \item{dev}{node deviance}
#' \item{loss}{node loss}
#' \item{split}{the split that derived each node}
#' \item{pclass}{the predicted classification for each node}
#' 
#' \code{print.grab} Summarises the node selected into the following:
#' Node ID, Number of observations, Number of sets, Number of predators,
#' Number of prey, deviance, expected loss and predicted class.
#' 
#' The output is printed to a .csv file
#' 
#' @seealso   \code{\link{print}}, \code{\link{diversity}}
#' 
#' @examples 
#' 
#' # Example 1: Printing output from dpart
#' # Load data
#' data(yftdiet)  
#' # Load the prey taxa data
#' data(PreyTaxonSort)
#' 
#' # Assigning prey colours for default palette
#' val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = TRUE)
#' node.colsY <- val$cols
#' dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes
#' 
#' # Fitting the classification tree
#' yft.dp <- dpart(Group ~ Lat + Lon + Year + Quarter + SST  + Length, 
#'                   data = dietPP, weights = W, minsplit = 10,
#'                                     cp = 0.001)
#'                                     plot(yft.dp, node.cols = node.colsY)
#'                                     summary(yft.dp)
#'                                     print(yft.dp, setID = "TripSetNo")
#'                                     
#' # Example 2: Printing nodes from the tree
#' # Load data
#' data(yftdiet)  
#' # Load the prey taxa data
#' data(PreyTaxonSort)
#' 
#' # Assigning prey colours for default palette
#' val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = TRUE)
#' node.colsY <- val$cols
#' dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes
#' 
#' # Fitting the classification tree
#' yft.dp <- dpart(Group ~ Lat + Lon + Year + Quarter + SST  + Length, 
#'                 data = dietPP, weights = W, minsplit = 10,
#'                                 cp = 0.001)
#'                                  
#' # Pruning the tree
#' yft.pr <- prune(yft.dp, se = 1)                   
#' # Exploring Nodes: This suite of graphics is interactive and therefore will not be run. 
#' # When run, the code will ask you to select a node for
#' # viewing.
#' 
#' # Exploring nodes of the tree - single page
#' \dontrun{
#' val <- grab(object = yft.pr, LatID = "Lat", LonID = "Lon", setID = "TripSetNo", 
#' #            node.cols = node.colsY, cex = 1, mapxlim = c(-125, -75), mapylim = c(0, 30),
#' #            mapcol = "gold3", pos = "topleft")
#' # val
#' }


print <- function(x, ...){
  UseMethod("print")
  
  
}

#' @rdname print
print.dpart <- function(x, setID = NULL, digits=getOption("digits"), file = "diet_tree_summary.csv", ...){
  
  object <- x
  ff <- object$frame
  node <- as.integer(row.names(ff))
  nobs <- ff$n
  npredators <- ff$wt
  dev <- ff$dev
  expLoss <- ff$dev/npredators
  
  val <- rpartco.dpart(object)
  temp <- rpart.branch(val$x, val$y, node, branch=1)
  nodeLR <- sort(c(1, temp$nodeL, temp$nodeR))
  nprey <- nsets <- pclass <- splits <- NULL
  for(i in 1:length(nodeLR)){
    subtree <- snip.dpart(object, nodeLR[i])
    nodesub <- as.numeric(row.names(subtree$frame))
    val <- rpartco.dpart(subtree)
    rn <- row.names(subtree$frame)[subtree$frame$var == "<leaf>"]
    pred.where <- names(table(subtree$where))[(1:length(rn))[rn == nodeLR[i]]]
    subpred <- predict(subtree, type = "prob", plot = FALSE)
    subpred.where <- subtree$where
    pred.node <- subpred[subpred.where == paste(pred.where),]
    pred.node.m <- apply(pred.node, 2, mean)
    nprey[i] <- length(pred.node.m[pred.node.m > 0])
    if(is.null(setID))
      nsets[i] <- NA
    else
      nsets[i]  <- length(unique(subtree$data[subpred.where ==
                                                as.integer(pred.where),][,setID]))
  }
  # Extracting the predicted class
  pclass <- levels(object$data$Group)[ff$yval]
  
  # Extracting the split information
  splits <- ff$var[order(node)]
  options(warn = -1) # suppresses an unnecessary warning message
  
  xc <- rpconvert(object)
  # xc <- object
  options(warn = 0)
  xc.sp <- xc$splits
  
  cuts <- vector(mode='character', length=nrow(xc$splits))
  temp <- xc$splits[ ,2L]
  for (i in 1L:length(cuts)) {
    if(temp[i] !=0){ 
      if (temp[i] == -1L)
        cuts[i] <-paste("<", format(signif(xc$splits[i,4], digits=digits)))
      else if (temp[i] ==1L)
        cuts[i] <-paste(">", format(signif(xc$splits[i,4], digits=digits)))
      else cuts[i]<- paste("splits as ",
                           paste(c("L", "-", "R")[xc$csplit[xc$splits[i,4], 1L:temp[i]]],
                                 collapse='', sep=''), collapse='')
    }
  }
  
  
  splitval <- vector("numeric", length = nrow(ff))
  for(i in 1:nrow(ff)){    
    m <- (row.names(xc$splits) == ff$var[i]) & 
      (xc$splits[,1] == ff$n[i]) 
    if(length(cuts[m]) !=0){
      splitval[i] <- cuts[m][1]
    }
    else
      splitval[i] <- NA
  }
  splitval <- splitval[order(node)]
  splitcomb <- paste(splits, splitval)
  splitcomb[is.na(splitval)] <- "<leaf>"
  
  sumdata <- data.frame(node = sort(node), nobs = nobs[order(node)], nsets = nsets,
                        npredators = npredators[order(node)], nprey = nprey, dev = dev[order(node)], 
                        loss = expLoss[order(node)], split = splitcomb, pclass = pclass[order(node)])
  
  write.csv(sumdata, file, row.names = FALSE)
  cat(paste("Written summary table out to file: ", file, "\n", sep = ""))
  
  invisible(sumdata)
  
}


#' @rdname print
print.grab <- function (x, ...) 
{
  
  object <- x
  cat("Node \t", "N \t", "No. of sets \t", "No. of predators \t", 
      "No. Prey \t", "Deviance \t", "Expected Loss \t Predicted Class\n")
  cat("------------------------------------------------------------------------------------------------------------------------------------\n")
  for (i in 1:nrow(object$nodeS))
    with(object$nodeS[i, ], cat(node,
                                "\t", nobs, "\t", nsets, "\t\t", npredators, "\t\t\t", 
                                nprey, "\t\t", round(dev, 2), "\t\t", round(loss, 4), "\t", as.character(pclass),
                                "\n"))
  
  # Printing to file
  write.csv(object$nodeS, "grab_node_res.csv", row.names = FALSE)
  cat("Written grabbed node summary table out to file: grab_node_res.csv\n")
  
  # Plotting bootstrap predictions
  if(any(names(object) == "nodestats")){
    nodestats <- object$nodestats
    options(warn = -1)
    par.old <- par()$mar
    par(mar = c(7, 4, 4, 2) + 0.1)
    plot(1:nrow(nodestats), nodestats$m, type = "n", axes = FALSE,
         ylim = c(0,1), xlab = "", ylab = "Bootstrapped Proportion")
    points(1:nrow(nodestats), nodestats$m, pch = 16)
    arrows(1:nrow(nodestats), nodestats$lci95, 1:nrow(nodestats), nodestats$uci95,
           angle = 90, code = 3, length = 0.1)
    box()
    axis(side = 2)
    mtext(side = 1, text = row.names(nodestats), adj = 1, line = 1,
          at = 1:nrow(nodestats), las = 2)
    
    title(main = paste("Node", object$nodeS$node))
    par(mar = par.old)
    options(warn = 0)
  }
  invisible()
}



