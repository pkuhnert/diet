#' Smooth Plots
#' 
#' @description Smooth plots used in \code{plot.diet}
#' 
#' @param x x
#' @param i i
#' @param SmXvar SmXvar
#' @param SmXdat SmXdat
#' @param LonID LonID
#' @param LatID LatID
#' @param projection projection
#' @param palette palette
#' @param too.far too.far
#' 
#' @importFrom "mgcv" "gam"
#' @import "sp"
#' @import "spaMM" 
#' @importFrom "raster" "raster"
#' @import "lattice"
#' @importFrom "latticeExtra" "layer"
#' @importFrom "rasterVis" "levelplot"
#' @importFrom "utils" "data"
#' @export
SmPlots <- function(x, i, SmXvar, SmXdat, LonID, LatID, projection, palette, too.far = 0.05){
  
  tmpdat <- data.frame(SmXdat[,i], x[,LonID], x[,LatID])
  names(tmpdat) <- c("y", "Longitude", "Latitude")
  
  res <- gam(y ~ s(Longitude, Latitude), data = tmpdat)
  
  # data(worldcountries)
  plotfit <- Vis.Gam(res, too.far = too.far, plot.type = "contour", add = FALSE)
  mat <- plotfit$mat
  names(mat) <- c("Longitude", "Latitude", "z")
  
  coordinates(mat) <- ~Longitude+Latitude  # coordinates are being set for the raster
  proj4string(mat) <- CRS(projection)  # projection is being set for the raster
  gridded(mat) <- TRUE  # a gridded structure is being set for the raster
  mat.raster <- raster(mat)  # the raster is being created
  
  zext <- range(mat$z, na.rm = TRUE)
  # Set up of plot
  
  p <- levelplot(mat.raster, maxpixels=4e6, margin=FALSE, cuts=length(palette)-1, col.regions=palette,
                 xlab = "", ylab = "", at = seq(zext[1], zext[2], length = 50),
                 main = SmXvar[i])

  # country layer
  country.layer <- layer(
    sp.polygons(worldcountries, fill="darkgray", col = "lightgray"),
    data=list(sp.polygons=sp.polygons, worldcountries=worldcountries) 
  )
  
  # points layer
  points.layer <- layer(
    panel.points(Longitude, Latitude, pch = 16, col = "black", cex = 0.6),
    data = tmpdat
  )
  

  smplot <- p + country.layer + points.layer

  list(res = res, smplot = smplot)
}