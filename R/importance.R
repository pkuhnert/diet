#' Variable Importance Ranking
#' 
#' @description Produces a variable importance ranking from a fitted \code{dpart} object.
#' 
#' @param x object of class \code{dpart}
#' 
#' @usage \method{importance}{dpart}(x)
#' 
#' @details A barplot of the variable importance ranking sorted from highest to lowest is produced 
#' in addition to the printed output.
#' 
#' @details Variable importance ranking for the variables used in the classification tree.
#' 
#' @references  Breiman, L., Friedman, J.H., Olshen, R.A. and Stone, C.J. (1984) Classification 
#' and Regression Trees. Wadsworth International.
#' 
#' @seealso \code{\link{dpart}}
#' 
#' @export
#' 
#' @examples 
#' 
#' # Assigning prey colours for default palette
#' val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = TRUE)
#' node.colsY <- val$cols
#' dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes
#' 
#' # Fitting the classification tree
#' yft.dp <- dpart(Group ~ Lat + Lon + Year + Quarter + SST  + Length, 
#'                  data = dietPP, weights = W, minsplit = 10,
#'                                     cp = 0.001)
#' yft.pr <- prune(yft.dp, se = 1)
#' vi <- importance(yft.pr)
importance <- function(x) UseMethod("importance")

#' @rdname importance
#' @export
importance.dpart <- function(x){

options(warn = -1)
sp <- x$splits
rn <- row.names(sp)
sp <- data.frame(sp)
sp$Variable <- rn
frame <- data.frame(x$frame)
nmsf <- names(frame)
frame <- frame[frame$var != "<leaf>", nmsf[-length(nmsf)]]

# Extract competing and surrogate splits
comp.split <- sp[sp$adj == 0,] # had changed this from subset(sp, adj == 0)
nc <- frame$ncompete
ns <- frame$nsurrogate
idfr <- 1:nrow(frame)
nc.id <- ns.id <- NULL
for(i in 1:nrow(frame)){
  nc.id <- c(nc.id, rep(idfr[i], nc[i]+1))
  ns.id <- c(ns.id, rep(idfr[i], ns[i]))
}
comp.split[,"ID"] <- nc.id
surr.split <- sp[sp$adj != 0,] ###subset(sp, adj != 0)
surr.split[,"ID"] <- ns.id
vi <- x$ordered
vi[vi == FALSE] <- 0
nms.vi <- names(vi)
for(i in 1:nrow(frame)){
  sub.comp <- comp.split[comp.split$ID == i,]  ##subset(comp.split, ID == i)
  sub.split <- surr.split[surr.split$ID == i,] ## subset(surr.split, ID == i)
  vi[sub.comp$Variable[1]] <- vi[sub.comp$Variable[1]] + sub.comp$improve[1]
  surr.adj <- sub.split$adj[match(sub.comp$Variable[-1], sub.split$Variable)]
  surr.adj[is.na(surr.adj)] <- 0
  vi[sub.comp$Variable[-1]] <- vi[sub.comp$Variable[-1]] +  sub.comp$improve[-1] *  surr.adj
  
}

varimp <- vi/max(vi)
varimp <- sort(varimp, decreasing = TRUE)

barplot(varimp, main = "Variable Importance")
print(varimp, digits = 2)

}



