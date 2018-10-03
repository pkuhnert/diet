#' Plot of the diet classification tree
#' 
#' @description 
#' The function produces a plot of the classification tree with terminal nodes
#' coloured according to the node colours provided. If \code{node.cols} is NULL, the standard
#'  classification tree is produced.
#'  
#'  @aliases plot
#'  
#' @param x object of class dpart.
#' @param y NULL (not used)
#' @param node.cols vector of colours with prey labels. Produced from the \code{apc} function.
#' @param keep.margins logical. Whether to keep the same plotting margins (TRUE) or reset 
#' them back to the default (FALSE). (Default: FALSE)
#' @param ... arguments to be passed to or from other methods
#'  
#' @details  The \code{keep.margins} should be left at the default. This argument is used
#' by the \code{grab} function to make the plot visually appealing.
#'  
#'  
#' @seealso  \code{\link{dpart}}, \code{\link{apc}}, \code{\link{rpart}}, \code{\link{plot.rpart}} 
#'  
#' @references Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
#' and Stable Isotope Data, Journal of Statistical Software, In Prep.
#'  
#' @examples 
#' # Load data
#' #data(yftdiet)  
#'  
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
#' #                    data = dietPP, weights = W, minsplit = 10, cp = 0.001)
#' # plot(yft.dp, node.cols = node.colsY)
#'                    
#' @import rpart.plot                    
#' @export
plot.dpart <- function(x, y = NULL, node.cols = NULL, keep.margins = FALSE, ...){


  # Plot 2: tree
  val <- prp(x, extra = 8, box.palette = as.list(node.cols), 
      shadow.col = "gray", split.prefix = "is ", 
      split.suffix = "?", nn = TRUE, split.cex = 1.5, cex = 0.6, yesno = TRUE,
      type = 2, branch.lwd = 3, legend.x = NULL, legend.y = NULL, roundint = FALSE, ...)
  
 


  invisible(val)
}




