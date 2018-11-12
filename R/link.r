#' Links Bagged Predictions to Tree
#' 
#' @description This function links the bagged predictions to a tree object
#' 
#' @param x object of class \code{bag}.
#' @param object tree object of class \code{dpart} that is used to link the bagged predictions to.
#' @param LatID column name for latitude.
#' @param LonID column name for longitude.
#' @param plot logical. Should a plot be produced.
#' @param oob logical. Should out of bag (oob) predictions be used to map back to the terminal nodes
#' of the tree.
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
#' 
#' # Assigning prey colours for default palette
#' val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = TRUE)
#' node.colsY <- val$cols
#' dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes
#' 
#' # Fitting the classification tree
#' yft.dp <- dpart(Group ~ Lat + Lon + Year + Quarter + SST  + Length, 
#'                    data = dietPP, weights = W, minsplit = 10,
#'                                      cp = 0.001)
#' yft.pr <- prune(yft.dp, se = 1)
#' plot(yft.pr, node.cols = node.colsY)
#' 
#' # Bagging
#' # Bagging with NO spatial bootstrapping
#' #yft.bag <- bagging(Group ~ Lat + Lon + Year + Quarter + SST  + Length,
#' #                       data = dietPP, weights = W, minsplit = 50,
#' #                       cp = 0.001, nBaggs = 500, predID = "TripSetPredNo")
#'                      
#' # Link the predictions back to terminal nodes of a tree and plot 
#' #ex.bag.l <- link(x = yft.bag, object = yft.pr, LatID = "Lat", LonID = "Lon",
#' #                   mapxlim = c(-125, -75), mapylim = c(0, 30), plot = TRUE, oob = TRUE)
#'             
#'  
#' @export                  
link <- function(x, object, LatID, LonID, plot = TRUE, oob = FALSE, orderN = FALSE) 
  UseMethod("link")

#' @rdname link
#' @importFrom "reshape2" "melt"
#' @import "abind"  
#' @importFrom "gridExtra" "marrangeGrob"
#' @export
link.bag <- function(x, object, LatID, LonID, 
                     plot = TRUE,  oob = FALSE, orderN = FALSE){
  
  # links the predictions to the terminal nodes of a tree model
  options(warn = -1)
  if (!inherits(x, "bag")) 
    stop("Not a bagged object")
  
  dat <- object$data
  nBaggs <- length(x$baggs)
  nodenms <- row.names(object$frame)[object$frame$var == "<leaf>"]
  
  
  bpred <- list()
  for(i in 1:nBaggs){
    bpred[[i]] <- data.frame(matrix(0, nrow = length(names(table(object$where))),
                                    ncol = ncol(x$pred[[i]])))
    row.names(bpred[[i]]) <- names(table(object$where))
    names(bpred[[i]]) <- names(x$pred[[i]])
  }
  options(warn = -1)
  for(i in 1:nBaggs){
    if(oob){  
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
  
  row.names(bpred.m) <- nodenms
  row.names(bpred.v) <- nodenms
  row.names(bpred.lci) <- nodenms
  row.names(bpred.uci) <- nodenms
  
  
  if(plot){

    bpred <- melt(data.frame(bpred.m))
    names(bpred)[2] <- "mean"
    bpred.lci <- melt(data.frame(bpred.lci))
    bpred.uci <- melt(data.frame(bpred.uci))
    # add a LCI column to bpred
    id <- match(row.names(bpred.lci), row.names(bpred))
    bpred$lci <- bpred.lci$value[id]
    # add a UCI column to bpred
    id <- match(row.names(bpred.uci), row.names(bpred))
    bpred$uci <- bpred.uci$value[id]
    # add node information
    bpred$node <- rep(row.names(bpred.m), ncol(bpred.m))
    
    bp <- list()
    
    un_node <- unique(bpred$node)
    for(i in 1:length(un_node)){
      tmp <- subset(bpred, node == un_node[i])
      bp [[i]] <- ggplot(tmp) + 
        geom_segment(aes(x = variable, y = lci, xend = variable, yend = uci),
                     size = 1, lineend = "butt") + geom_point(aes(x = variable, y = mean), col = "white", size = 1, shape = "-") +
                     xlab("") + ylab("Proportion") + theme_bw() + ggtitle(paste("Node ", un_node[i])) +
        theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) 
    }


    m_bp <- marrangeGrob(grobs = bp, nrow = 3, ncol = 3, top = "Node Predictions with 95% Bootstrapped CIs")
    # print to screen
    m_bp
    # save to file
    ggsave("Node_Predictions.pdf", m_bp, width = 8, height = 6)
    cat("Figure saved to file: Node_Predictions.pdf")
    val <- list(m = bpred.m, v = bpred.v, lci95 = bpred.lci, uci95 = bpred.uci, m_bp = m_bp)
    
  }
  else
    val <- list(m = bpred.m, v = bpred.v, lci95 = bpred.lci, uci95 = bpred.uci, m_bp = NULL)
  
  class(val) <- c("link")
  
  val
  
} 



