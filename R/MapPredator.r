MapPredator <- function(x, LonID, LatID, mapxlim, mapylim, database, PredSpID){
  
  expl2 <- mappoints.data(x[,c(LonID,LatID)], xlim = mapxlim, ylim = mapylim,  
                          database = database)
  XYdata <- x[,c(LonID, LatID)]
  Longitude <- XYdata[,1]
  Latitude <- XYdata[,2]
  df <- data.frame(Longitude, Latitude, Predator = x[,PredSpID])
  expl3 <- expl2 + geom_point(data = df, aes(Longitude, Latitude, colour = Predator)) +
    ggtitle("Distribution of Samples by Predator")
  
  expl3
}