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





