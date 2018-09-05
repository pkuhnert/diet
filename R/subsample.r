# Takes a random sub-sample, conditional on ID from a dataset
subsample <- function(dat, ID, n){

  x <- NULL
  unID <- unique(ID)
  t.ID <- table(ID)
  if(all(t.ID <= n))
     x <- dat
  else{
    # do subsampling
    for(i in 1:length(unID)){
      subdat <- subset(dat, ID == unID[i])
      if(is.null(n))
         x <- rbind(x, subdat[sample(1:nrow(subdat), size = nrow(subdat), replace = TRUE),])
      else
         if(nrow(subdat) <= n)
            x <- rbind(x, subdat[sample(1:nrow(subdat), size = nrow(subdat), replace = TRUE),])
      else
         x <- rbind(x, subdat[sample(1:nrow(subdat), size = n, replace = TRUE),])
  }
  }
  

  
  x

}


