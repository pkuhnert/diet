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
       plot(v.res, pch = 16, max.dist = max(v.res$u)/2)
    }

    resids
        
        
}
