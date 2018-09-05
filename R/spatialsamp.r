# Implementing a spatial bootstrap

spatialsamp <- function(x, LonID, LatID, sizeofgrid = 5, ID = NULL, nsub = NULL, Plot = FALSE){


   lon.seq <- seq(min(x[,LonID], na.rm = TRUE), max(x[,LonID], na.rm = TRUE),
      by = sizeofgrid)
   lat.seq <- seq(from = min(x[,LatID], na.rm = TRUE), to = max(x[,LatID], na.rm = TRUE),
      by = sizeofgrid)

   if(Plot){
      plot(x[,LonID], x[,LatID], pch = 16, cex = 0.7)
      abline(h = lat.seq, lty = 3, col = "grey")
      abline(v = lon.seq, lty = 3, col = "grey")
    }

   gx <- cut(x[,LonID], lon.seq, include.lowest = TRUE)
   gy <- cut(x[,LatID], lat.seq, include.lowest = TRUE)
   ugx <- unique(gx)
   ugy <- unique(gy)
   subsamp <- NULL
   k <- 1

   for(i in 1:length(ugx)){
       for(j in 1:length(ugy)){

           ids <- (1:nrow(x))[gx %in% ugx[i] & gy %in% ugy[j]]
           samp <- x[ids,]
           if(Plot)
              with(samp, points(Lon, Lat, pch = 16, col = k, cex = 0.7))

           if(is.null(nsub)){
              Ssamp <- sample(1:nrow(samp), nrow(samp), replace = TRUE)
              subsamp <- rbind(subsamp, samp[Ssamp,])
           }
           else{
              Ssamp <- subsample(dat = samp, ID = samp[,ID], n = nsub)
              subsamp <- rbind(subsamp, Ssamp)
           }
           k <- k+1
       }
    }

   subsamp

}
        
        
