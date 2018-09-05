read.pp <- function (object = NULL, filenm = NULL, 
                     labels = list(PredatorL = "PredID",TripSetL = "TripSetNo", SpeciesL = "Species", 
                                   FullnessL = "Fullness",DateL = "Date", WeightL = "Weight", 
                                   PreyGrpL = "PreyGrp"), 
                     Datef = "%d/%m/%Y", p = 0.01, Xvars, sp2omit = NULL) 
{
  if (is.null(object) & is.null(filenm)) 
    stop("Either a filename or object needs to be provided.")
  if (missing(labels)) 
    stop("Labels specifying the column names for extract the diet indices have not\n          been supplied.")

  
  if (is.null(object)) 
    x <- read.csv(filenm)
  else x <- object
  PredatorL <- labels$PredatorL
  TripSetL <- labels$TripSetL
  SpeciesL <- labels$SpeciesL
  FullnessL <- labels$FullnessL
  WeightL <- labels$WeightL
  PreyGrpL <- labels$PreyGrp
  DateL <- labels$DateL
  dat <- x[order(x[, PredatorL], x[, SpeciesL]), ]
  dat[, PreyGrpL] <- as.factor(as.vector(dat[, PreyGrpL]))
  d <- dat[, DateL]
  d <- as.POSIXct(as.vector(d), format = Datef, tz = "GMT")
  dat[, "Date"] <- d
  f <- dat[, FullnessL]
  dat[, "Fullness"] <- f
 
  Xmat <- data.frame(matrix(NA, nrow = length(unique(dat[, 
                                                         PredatorL])), ncol = length(Xvars) + 4))
  names(Xmat) <- c(PredatorL, TripSetL, DateL, FullnessL, Xvars)
  Ymat <- data.frame(matrix(NA, nrow = length(unique(dat[, 
                                                         PredatorL])), ncol = length(unique(dat[, PreyGrpL]))))
  names(Ymat) <- levels(dat[, PreyGrpL])
  uPid <- unique(dat[, PredatorL])
  row.names(Ymat) <- uPid
  Ymat[is.na(Ymat)] <- 0
  for (i in 1:length(uPid)) {
    tmp <- dat[dat[, PredatorL] == uPid[i], ]
    m <- match(as.vector(tmp[, PreyGrpL]), names(Ymat))
    wt <- tapply(tmp[, WeightL], m, sum)
    Ymat[i, as.integer(names(wt))] <- wt
    Xmat[i, ] <- tmp[, match(c(PredatorL, TripSetL, DateL, 
                               FullnessL, Xvars), names(tmp))][1, ]
  }
  # Keeping factor levels for variables
  Xcheck <- dat[,Xvars]
  for(i in 1:length(Xvars)){ 
    Xvar.logical <- is.factor(Xcheck[,i])
    if(Xvar.logical){
      Xmat[,match(names(Xcheck)[i], names(Xmat))] <- factor(Xmat[,match(names(Xcheck)[i], names(Xmat))])
      levels(Xmat[,match(names(Xcheck)[i], names(Xmat))]) <- levels(dat[,Xvars][,i])
    }
  }
  
  Xmat$Date <- dat$Date[match(Xmat[, PredatorL], dat[, PredatorL])]
  diet.matrix <- cbind(Xmat, Ymat)
  diet.pp <- read.dm(object = diet.matrix, 
                     labels = list(DateL = "Date",FullnessL = "Fullness"), 
                     diet.ind.start = ncol(Xmat) + 1, p = p, sp2omit = sp2omit)

  
 diet.pp
}

