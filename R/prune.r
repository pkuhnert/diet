#' Cost-complexity pruning of a dpart object
#' 
#' @description This function prunes a fitted \code{dpart} object based on either the standard
#' error or number of splits required.
#' 
#' @param tree a legitimate tree object of class \code{dpart}.
#' @param se numeric. A standard error used to prune. 
#' @param nsplits  numeric. Number of splits to prune to.
#' @param \dots arguments to be passed to or from other methods
#' 
#' @usage   \method{prune}{dpart}(tree, se, nsplits, ...)
#' 
#' @details   A new \code{dpart} object that is pruned based on the number of standard errors
#' or number of splits.
#' 
#' @references   Breiman, L., Friedman, J.H., Olshen, R.A. and Stone, C.J. (1984) Classification 
#' and Regression Trees. Wadsworth International.
#' 
#' @seealso \code{\link{dpart}}
#' 
#' @export
#' 
#' @examples 
#' 
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
#' #                   data = dietPP, weights = W, minsplit = 10,
#' #                                     cp = 0.001)
#' #yft.pr <- prune(yft.dp, se = 1)
#' #plot(yft.pr, node.cols = node.colsY)


prune <- function(tree, ...){
      UseMethod("prune")
}

#' @rdname prune
prune.dpart <- function(tree, se, nsplits, ...){
  
  
  tree.cp <- select.tree(tree, se = se, nsplits = nsplits)
  prune.rpart(tree, cp = tree.cp, ...)
  
}

