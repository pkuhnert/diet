#' MapPreyPredator
#' 
#' @description Internal function used in \code{plot.diet}
#' 
#' @param x x
#' @param PredSpID PredSpID
#' @param LonID LonID
#' @param LatID LatID
#' @param mapxlim mapxlim
#' @param mapylim mapylim
#' @param database database
#' 
#' @importFrom "ggplot2" "ggplot" "ggplot_gtable" "aes_string" "geom_point" "geom_bar"
#' @importFrom "ggplot2" "geom_histogram"

MapPreyPredator <- function(x, PredSpID, LonID, LatID, mapxlim, mapylim,
                database){
  
  
  unpID <- unique(x[,PredSpID])
  if(length(unpID) != 1){
    for(i in 1:length(unpID)){
      # Extract predator and plot for each prey
      pdat <- x[x[,PredSpID] == unpID[i],]
      unpreyID <- levels(pdat$Group)
      map2 <- list()
      for(j in 1:length(unpreyID)){
        prdat <- pdat[pdat$Group == unpreyID[j],]
        
        map2[[j]] <- mappoints(prdat[,c(LonID,LatID)], xlim = mapxlim, ylim = mapylim,  
                               database = database, 
                               gtitle = paste("Predator: ", unpID[i], "\nPrey: ", unpreyID[j], sep = ""))
      }
    }
    expl5 <- do.call(grid.arrange, map2)
    
  }
  else
    expl5 <- NULL
  expl5
}