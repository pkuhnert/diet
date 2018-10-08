#' Residual Analysis
#' 
#' @description Examines the residuals from the fit of the model to determine whether spatial
#' dependence exists. Operates on objects of class \code{dpart}.
#' 
#' @param object object of class \code{dpart}.
#' @param LonID column name for longitude.
#' @param LatID column name for latitude.
#' @param predID predator ID
#' @param plot logical. Whether plotting should be performed
#' 
#' @usage 
#' \method{resids}{dpart}(object, LonID, LatID, predID, plot = TRUE)
#' \method{resids}{bag}(object, LonID, LatID, predID, plot = TRUE)
#' 
#' 
#' @details   Produces a variogram of the residuals to check for spatial dependence.
#' 
#' @return  Outputs the residuals and a variogram plot of the residuals for checking.
#' 
#' @seealso objects to See Also as \code{\link{variog}}
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
#' #                   data = dietPP, weights = W, minsplit = 10,
#' #                                     cp = 0.001)
#'                                     
#' # Pruning the tree
#' #yft.pr <- prune(yft.dp, se = 1)                   
#' 
#' # Checking residuals
#' #yft.resid <- resids(yft.pr, LonID = "Lon", LatID = "Lat", 
#' #     predID = "TripSetPredNo", plot = TRUE) # need to compute resids from bootstrapping
#'     
#' # Bagging
#' # Bagging with NO spatial bootstrapping (nBaggs set to something small)
#' #yft.bag <- bagging(Group ~ Lat + Lon + Year + Quarter + SST  + Length,
#' #                     data = dietPP, weights = W, minsplit = 50,
#' #                     cp = 0.001, nBaggs = 10, predID = "TripSetPredNo")
#'                     
#' # Checking for spatial dependence
#' #yft.bag.resid <- resids(yft.bag, LonID = "Lon", LatID = "Lat", 
#' #predID = "TripSetPredNo", plot = TRUE) # need to compute resids from bootstrapping
#'           
#'                     
#' @import "geoR"
#'                                         
#' @export                     


resids <- function(object, LonID, LatID, predID, plot = TRUE) 
  UseMethod("resids")

#' @rdname resids
#' @import "geoR"
#' @export
resids.dpart <- function(object, LonID, LatID, predID, plot = TRUE){
  
  
  if (!inherits(object, "dpart")) 
    stop("Not diet object")
  
  Omat <- formOmat(object$data, ID = predID)
  pp <- predict(object, type = "prob", pred.type = "predator",
                predatorID = predID, plot = FALSE)
  resids <- Distance(O = Omat[,(ncol(Omat)-ncol(pp)+1):ncol(Omat)], 
                     P = pp, type = "Hellinger")
  
  
  
  browser()
  if(plot){
    spdat <- data.frame(Lon = Omat$Lon, Lat = Omat$Lat)
    id <- (1:nrow(spdat))[apply(spdat, 1, function(x) any(is.na(x)))]
    if(length(id) == 0)
      v.res <- variog(coords = spdat, data = resids$Dist)
    else
      v.res <- variog(coords = spdat[-id,], data = resids$Dist[-id])
 
    v_df <- data.frame(distance = v.res$u, semivariance = v.res$v)
    semi_plot <- ggplot(v_df, aes(distance, semivariance)) + geom_point(size = 2) + xlim(0, max(v_df$distance)/2) +
      theme_bw() + ggtitle("Variogram of Residuals")
    print(semi_plot)
  }
  
  resids
  
  
}

#' @rdname resids
#' @import "geoR"
#' @export
resids.bag <- function(object, LonID, LatID, predID, plot = TRUE){
  
  
  if (!inherits(object, "bag")) 
    stop("Not an bagged object")
  
  pred <- do.call(abind, c(object$predP, along = 3))
  pred <- apply(pred, c(1,2), mean)
  
  Omat <- formOmat(object$data, ID = predID)
  
  resids <- Distance(O = Omat[,(ncol(Omat)-ncol(pred)+1):ncol(Omat)], 
                     P = pred, type = "Hellinger")
  
  
  
  
  if(plot){
    spdat <- data.frame(Lon = Omat[,LonID], Lat = Omat[,LatID], z = resids$Dist)
    spdat <- na.omit(spdat)
    v.res <- variog(coords = spdat[,c("Lon", "Lat")], data = spdat$z)
    v_df <- data.frame(distance = v.res$u, semivariance = v.res$v)
    semi_plot <- ggplot(v_df, aes(distance, semivariance)) + geom_point(size = 2) + xlim(0, max(v_df$distance)/2) +
      theme_bw() + ggtitle("Variogram of Residuals")
    print(semi_plot)
    
  }
  
  resids
  
  
}

