#' MapPrey
#' 
#' @import "ggplot2"
MapPrey <- function(x, LonID, LatID, mapxlim, mapylim, database){
  
  unpreyID <- levels(x$Group)
  map1 <- list()
  for(i in 1:length(unpreyID)){
     pdat <- x[x$Group == unpreyID[i],]
     map1[[i]] <- mappoints(pdat[,c(LonID,LatID)], xlim = mapxlim, ylim = mapylim,  
                         database = database, gtitle = unpreyID[i])
  }
  expl4 <- do.call(grid.arrange, map1)
  expl4
}