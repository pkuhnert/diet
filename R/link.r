#' Links Bagged Predictions to Tree
#' 
#' @description This function links the bagged predictions to a tree object
#' 
#' @param x object of class \code{bag}.
#' @param object tree object of class \code{dpart} that is used to link the bagged predictions to.
#' @param LatID column name for latitude.
#' @param LonID column name for longitude.
#' @param mapxlim vector of upper and lower limits for the x-axis of the map.
#' @param mapylim vector of upper and lower limits for the y-axis of the map.
#' @param plot logical. Should a plot be produced.
#' @param oob logical. Should out of bag (oob) predictions be used to map back to the terminal nodes
#' of the tree.
#' @param mfrow a vector of the form \code{c(nr,nc)} where \code{nr} equals the number of rows and
#' \code{nc} equals the number of columns for plotting. (defaults to c(2,2))
#' @param orderN  logical. Should the node outputs be ordered.
#' 
#' 
#' 
#' @return a list containing the following items:
#' \item{m}{matrix of mean proportions}
#' \item{v}{matrix of variances}
#' \item{lci95}{matrix of lower 95\% confidence intervals}
#' \item{uci95}{matrix of upper 95\% confidence intervals}
#' 
#' @details   The direction of a zero-length arrow is indeterminate, and hence so 
#' is the direction of the arrowheads. To allow for rounding error, 
#' arrowheads are omitted (with a warning) on any arrow of length less 
#' than 1/1000 inch. See \code{arrows} for more details.
#' 
#' @references   Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
#' and Stable Isotope Data, Journal of Statistical Software, In Prep. 
#' 
#' Kuhnert PM, Kinsey-Henderson A, Bartley R, Herr A (2010)
#' Incorporating uncertainty in gully erosion calculations using
#' the random forests modelling approach. Environmetrics
#' 21:493-509. doi:10.1002/env.999
#' 
#' @seealso   \code{\link{bagging}}
#' @examples 
#' # Load data
#' #data(yftdiet)  
#' # Load the prey taxa data
#' #data(PreyTaxonSort)
#' 
#' # Assigning prey colours for default palette
#' #val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = TRUE)
#' #node.colsY <- val$cols
#' #dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes
#' 
#' # Fitting the classification tree
#' #yft.dp <- dpart(Group ~ Lat + Lon + Year + Quarter + SST  + Length, 
#' #                   data = dietPP, weights = W, minsplit = 10,
#' #                                     cp = 0.001)
#' #yft.pr <- prune(yft.dp, se = 1)
#' #plot(yft.pr, node.cols = node.colsY)
#' 
#' # Bagging
#' # Bagging with NO spatial bootstrapping
#' #yft.bag <- bagging(Group ~ Lat + Lon + Year + Quarter + SST  + Length,
#' #                      data = dietPP, weights = W, minsplit = 50,
#' #                      cp = 0.001, nBaggs = 500, predID = "TripSetPredNo")
#'                      
#' # Link the predictions back to terminal nodes of a tree and plot 
#' #ex.bag.l <- link(x = yft.bag, object = yft.pr, LatID = "Lat", LonID = "Lon",
#' #                  mapxlim = c(-125, -75), mapylim = c(0, 30), plot = TRUE, oob = TRUE, 
#' #                  mfrow = c(2,2))
#'                    
#' @export                  
link <- function(x, object, LatID, LonID, mapxlim, mapylim,
   plot = TRUE, oob = FALSE, mfrow = c(2,2), orderN = FALSE) 
  UseMethod("link")

#' @rdname link
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
      #where <- pred.rpart(object, rpart.matrix(x$data[paste(x$oob[[i]]),]))
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
   # if(dev.cur() == 1) windows(8,8, record = TRUE)
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



