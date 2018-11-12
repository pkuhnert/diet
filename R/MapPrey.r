#' MapPrey
#' 
#' @description Internal function used by \code{plot.diet}
#' 
#' @param x x
#' @param LonID LonID
#' @param LatID LatID
#' @param mapxlim mapxlim
#' @param mapylim mapylim
#' @param database database
#' 
#' @importFrom "ggplot2" "ggplot" "ggplot_gtable" "aes_string" "geom_point" "geom_bar"
#' @importFrom "ggplot2" "geom_histogram"


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