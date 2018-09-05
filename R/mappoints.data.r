mappoints.data <- function(XY.data, xlim, ylim, database = 'world', gtitle = NULL){

  make_labels <- function(value, loc){
    
    paste(value, "º", loc, sep = "")
    
  }


  if(missing(XY.data))
     stop("Longitude and Latitude required when producing a map.\n")
  
  Longitude <- XY.data[,1]
  Latitude <- XY.data[,2]
  if(missing(xlim))
       xlim <- round(range(XY.data[,1], na.rm = TRUE))

  if(missing(ylim))
     ylim <- round(range(XY.data[,2], na.rm = TRUE))

  xseq <- seq(xlim[1], xlim[2], 15)
  xloc <- ifelse(xseq < 0, "W", "E")
  xloc[xseq ==0] <- ""
  
  
  yseq <- seq(ylim[1], ylim[2], 15)
  yloc <- ifelse(yseq < 0, "S", "N")
  yloc[yseq == 0] <- ""
  
  mapWorld <- borders(database, colour="gray50", fill="gray50") 
  mp <- ggplot() +  mapWorld  + xlab("Longitude") + 
    ylab("Latitude") + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
  mp2 <- mp + geom_point(aes(Longitude, Latitude), colour = "black", data = XY.data, 
                         size = 1.5)
  mp3 <- mp2 + scale_x_continuous(limits = xlim, breaks = xseq, 
                                  labels = make_labels(value = xseq, loc = xloc)) +
    scale_y_continuous(limits = ylim, breaks = yseq, 
                       labels = make_labels(value = yseq, loc = yloc))


  
  if(is.null(gtitle))
    mp4 <- mp3 + ggtitle("Samples")
  else
    mp4 <- mp3 + ggtitle(gtitle)
  mp4


}