#' 
#' Produce a Map of a Region with Data overlayed
#' 
#' @description   
#' Utilises the \code{map} function to produce a map of the study area based on
#' locations given.
#' 
#' @usage 
#' mappoints(XY.data, xlim, ylim, database = 'world', gtitle)
#'                            
#' @param XY.data matrix of longitude and latitude co-ordinates
#' @param xlim  map x limits
#' @param ylim map y limits
#' @param database either "world" or "world2" depending on the area being plotted. Defaults to "world"
#' @param gtitle Title of Plot             
#' 
#' @details 
#' Calls on the \code{map} function for plotting.
#' 
#' 
#' @references 
#' Richard A. Becker, and Allan R. Wilks, Maps in S, AT\&T Bell Laboratories 
#' Statistics Research Report [93.2], 1993. 
#' 
#' Richard A. Becker, and Allan R. Wilks, Constructing a Geographical Database, 
#' AT&T Bell Laboratories Statistics Research Report [95.2], 1995. 
#' 
#' @examples 
#' # Load dataset
#' data(yftdiet)
#' 
#' # Produce map of the data
#' mappoints(XY.data = yftdiet[,c("Lon", "Lat")])
#' #title(main = "YFT Data Distribution")
#' 
#' @import "maps"
#' @importFrom "ggplot2" "ggplot" "ggplot_gtable" "aes_string" "geom_point" "geom_bar"
#' @importFrom "ggplot2" "geom_histogram"

#' 
#' @export
mappoints <- function(XY.data, xlim, ylim, database = 'world', gtitle = NULL)
  UseMethod("mappoints")

#' @rdname mappoints
#' @export
mappoints.data <- function(XY.data, xlim, ylim, database = 'world', gtitle = NULL){
  
  make_labels <- function(value, loc){
    
   # paste(value, "º", loc, sep = "")
    paste(value, "\u00B0", loc, sep = "")
    
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


