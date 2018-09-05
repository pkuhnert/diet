grabmulti <-function(object, LatID, LonID, setID = NULL, node.cols = NULL, cex = 0.8,
         mapxlim = NULL, mapylim = NULL, mapcol = "gold3", n = nrow(object$frame), 
                     database = 'world', ...) UseMethod("grabmulti")