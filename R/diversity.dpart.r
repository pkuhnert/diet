diversity.dpart <- function(object, newdata = NULL, too.far = 0.1, LatID, LonID,
                            projection = "+proj=longlat +datum=WGS84", maxpixels = 1e5,
                            col = "lightgrey", fill = "lightgrey", zlim = c(0,1)){
 
   options(warn = -1)
   data(worldcountries)  # loading from package maps
   palette <- spaMM.colors()

   if(missing(LatID) | missing(LonID))
     stop("LatID and LonID column names are not specified.\n")
   
   if(is.null(newdata))
     newdata <- object$data
   
   ff <- object$frame
   npredators <- ff$wt
   dev <- ff$dev
   expLoss <- dev/npredators
   node <- as.integer(row.names(ff))
   
   
   D.matrix <- data.frame(node = node, D = expLoss) 
   
   val <- rpartco.dpart(object)
   temp <- rpart.branch(val$x, val$y, node, branch=1)
   x1 <- c(temp$x)[seq(1, length(c(temp$x)), by = 5)]
   y1 <- c(temp$y)[seq(1, length(c(temp$y)), by = 5)]
   x2 <- c(temp$x)[seq(4, length(c(temp$x)), by = 5)]
   y2 <- c(temp$y)[seq(4, length(c(temp$y)), by = 5)]
   x <- c(x1, x2)
   y <- c(y1, y2)
   nodeLR <-  as.integer(row.names(object$frame[object$frame$var == "<leaf>",]))
   d.bks <- seq(0, 1, length = 10)
   d.cols <- as.factor(d.bks)
   levels(d.cols) <- topo.colors(10)
   pdat <- NULL
   pred <- rpart:::pred.rpart(object, rpart.matrix(newdata))
   dd <- with(ff, dev/npredators)[pred]
   pdat <- data.frame(Lat = newdata[,LatID], Lon = newdata[,LonID], D = dd)
   names(pdat)[1:2] <- c("Latitude", "Longitude")
   pdat$Df <- cut(pdat$D, breaks = seq(0, 1, length = 10), include.lowest = TRUE)
   cols <- topo.colors(length(levels(pdat$Df)))
   pdat$cols <- cols[pdat$Df]
   row.names(pdat) <- row.names(newdata)

   

   # spatial prediction of diversity
   fit <- gam(D ~ s(Longitude, Latitude), data = pdat)
   plotfit <- Vis.Gam(fit, too.far = too.far, plot.type = "contour", 
                      add = FALSE, type = "response")
   mat <- plotfit$mat
   names(mat) <- c("Longitude", "Latitude", "Diversity")
   
   coordinates(mat) <- ~Longitude+Latitude  # coordinates are being set for the raster
   proj4string(mat) <- CRS(projection)  # projection is being set for the raster
   gridded(mat) <- TRUE  # a gridded structure is being set for the raster
   mat.raster <- raster(mat)  # the raster is being created
   
   
   # Set up of plot
   p <- levelplot(mat.raster, maxpixels=1e5, margin=FALSE, cuts=length(palette)-1, 
                  col.regions=palette, at = seq(0, 1, length.out = 20),
                  xlab = "", ylab = "")

   # country layer
   country.layer <- layer(
     sp.polygons(worldcountries, fill=fill, col = col),
     data=list(sp.polygons=sp.polygons, worldcountries=worldcountries, 
               fill=fill, col = col) 
   )
   
   # points layer
   points.layer <- layer(
     panel.points(Longitude, Latitude, pch = 16, col = "black", cex = 0.6),
     data = pdat
   )
  
   
   print(p + country.layer + points.layer)
}