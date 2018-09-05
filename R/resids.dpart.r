resids.dpart <- function(object, LonID, LatID, predID, plot = TRUE){

  
    if (!inherits(object, "dpart")) 
        stop("Not diet object")
  
    Omat <- formOmat(object$data, ID = predID)
    pp <- predict(object, type = "prob", pred.type = "predator",
                            predatorID = predID, plot = FALSE)
    resids <- Distance(O = Omat[,(ncol(Omat)-ncol(pp)+1):ncol(Omat)], 
                       P = pp, type = "Hellinger")

    
    
 
    if(plot){
       spdat <- data.frame(Lon = Omat$Lon, Lat = Omat$Lat)
       id <- (1:nrow(spdat))[apply(spdat, 1, function(x) any(is.na(x)))]
       if(length(id) == 0)
         v.res <- variog(coords = spdat, data = resids$Dist)
       else
         v.res <- variog(coords = spdat[-id,], data = resids$Dist[-id])
       plot(v.res, pch = 16, max.dist = max(v.res$u)/2)
    }
 
    resids
        
        
}
