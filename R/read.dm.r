read.dm <- function (object = NULL, filenm = NULL, labels = list(DateL = "Date", 
                    FullnessL = "Fullness"), Datef = "%d/%m/%Y", diet.ind.start, 
                     p = 0.01, sp2omit = NULL, predsp = NULL) 
{

  DateL <- labels$DateL
  FullnessL <- labels$FullnessL
  
  if (is.null(object) & is.null(filenm)) 
    stop("Either a filename or object needs to be provided.")
  
 
  if (is.null(object)) 
    dat <- read.csv(filenm)
  else dat <- object
  
  
  if(!is.null(predsp)){ 
     unpred <- as.vector(unique(dat[,predsp]))
     npred <- length(unpred)
     PreyData <- NULL
     for(i in 1:npred){
        cat("\nReading in ", unpred[i], " data. \n")

        tmp <- try(readpred(dat = dat[dat[,predsp] == unpred[i],], 
                        diet.ind.start = diet.ind.start, sp2omit = sp2omit, 
                        FullnessL = FullnessL, DateL = DateL, Datef = Datef, p = p), silent = TRUE)
        if(any(class(tmp) == "try-error"))
           cat("... no data for this species.\n")
        else
          PreyData <- rbind(PreyData, tmp)
        
     }
  }
  else
    PreyData <- readpred(dat = dat, diet.ind.start = diet.ind.start, 
                         sp2omit = sp2omit, FullnessL = FullnessL, DateL = DateL, Datef = Datef, p = p)
  
  
  oldClass(PreyData) <- c("diet", class(PreyData))
  
  writepn.csv(PreyData)
  
  PreyData
}