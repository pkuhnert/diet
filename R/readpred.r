#' Internal function used by \code{read.pp} and \code{read.dm}
#' 
#' @description Internal function that extracts the species information
#' 
#' @param dat A compositional data matrix containing covariate 
#' information and diet composition for each predator
#' @param diet.ind.start index where diet compositional data starts
#' @param sp2omit species to omit
#' @param FullnessL Fullness column name used to identify full stomachs for analysis
#' @param DateL Date column name
#' @param Datef date format
#' @param p Proportion of prey to omit. 


readpred <- function(dat, diet.ind.start, sp2omit, FullnessL, DateL, Datef, p){ 
 
dat <- dat[!apply(dat[, 1:(diet.ind.start - 1)], 1, function(x) all(is.na(x) | x == "")), ]
if (!is.null(sp2omit)) {
  if (any(is.na(match(sp2omit, names(dat))))) 
    warning("Species to omit is not listed. Therefore not removing any species from the list.")
  else {
    cat("Removing the following prey species ...\n")
    print(sp2omit)
    dat <- dat[, -match(sp2omit, names(dat))]
  }
}
dat <- dat[!apply(dat, 1, function(x) all(is.na(x))), ]
datWt <- dat[FullnessL != 0, ]
d <- datWt[, DateL]
if (!any(class(d) == "POSIXct")) {
  d <- as.POSIXct(as.vector(d), format = Datef, tz = "GMT")
  datWt$Date <- d
}
Ymat <- datWt[, diet.ind.start:ncol(dat)]
Ymat.nms <- names(Ymat)
Ymat[is.na(Ymat)] <- 0
Ymat.sum <- apply(Ymat, 1, function(x) sum(x))
datF <- datWt[Ymat.sum != 0, ]
datF <- datWt
Ymat <- datF[, diet.ind.start:ncol(dat)]
Ymat[is.na(Ymat)] <- 0
Ymat.sum <- apply(Ymat, 1, function(x) sum(x))
Ymat.p <- Ymat/Ymat.sum
Ymat.pvec <- apply(Ymat.p, 2, sum)/nrow(Ymat.p)

largeprey <- names(Ymat.pvec)[Ymat.pvec > p]
preyID <- match(largeprey, Ymat.nms)
if(length(preyID) == 1)
  Ymat <- data.frame(matrix(Ymat[, preyID], byrow = TRUE, nrow = nrow(Ymat)))
else
  Ymat <- Ymat[,preyID]
names(Ymat) <- largeprey
Ymat.sum <- apply(Ymat, 1, function(x) sum(x))
Ymat.p <- Ymat/Ymat.sum
datFprey <- cbind(datF[, 1:(diet.ind.start - 1)], Ymat.p)
Xvars <- datFprey[, 1:(diet.ind.start[1] - 1)]

PreyData <- NULL
for (i in 1:length(largeprey)) {
  PreyData <- rbind(PreyData, with(Ymat.p, cbind(Xvars, 
                                                 get(largeprey[i]), Group = largeprey[i])))
}
names(PreyData)[ncol(PreyData)-1] <- "W"
PreyData <- PreyData[PreyData$W > 0 & !is.na(PreyData$W),]

class(PreyData) <- c("data", class(PreyData))
for (i in 1:ncol(PreyData)) {
  if (is.factor(PreyData[, i])) 
    PreyData[, i] <- as.factor(as.vector(PreyData[, i]))
}
# Coerce Group column to factor
PreyData$Group <- as.factor(PreyData$Group)

 PreyData

}