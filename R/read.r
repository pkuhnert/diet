#' Read in Diet Matrix
#' 
#' @description Read in diet compositional data in one of two formats. Data 
#' stored as a compositional matrix can be read in using the \code{read.dm} 
#' function. Data stored in a predator-prey format can be read in using the 
#' \code{read.pp} function. 
#' 
#' @aliases read, read.dm
#' 
#' @param object matrix. A compositional data matrix containing covariate 
#' information and diet composition for each predator (optional).
#' @param filenm file name. A file name specifying the location of a comma 
#' delimited file that contains the covariate information and diet composition 
#' for each predator (optional).
#' @param labels a list containing the following: 
#' \code{DateL} Date column name: default label is "Date".
#' \code{FullnessL} Fullness column name used to identify full stomachs for analysis: 
#' default label is "Fullness".
#' \code{PredatorL} Predator column name: default label is "PredID".
#' \code{TripSetL} Sampling column name: default label is "TripSetNo".
#' \code{SpeciesL} Species column name: default label is "Species". This may 
#' have the same label as PreyGrpL.
#' \code{WeightL} Weight column name: default label is "Weight". This may 
#' represent weights, frequencies or binary data.
#' \code{PreyGrpL} Prey group column name: default label is "PreyGrp". This may 
#' have the same label as SpeciesL.
#' @param Datef date format. See \code{\link{strptime}} for format options.
#' @param diet.ind.start integer. Column number where the diet composition data begins.
#' @param p numeric. Proportion of prey to omit. If the composition is represented as a 
#' wet weight, then \code{p} represents the proportion of prey contributing to the wet 
#' weight consumption of a predator. If the composition is represented as a frequency or 
#' binary outcome, then \code{p} represents the proportion of prey found in a predator. (default: 0.01)
#' @param sp2omit a vector of prey names to omit from the analysis (default: NULL)
#' @param predsp Predator species names to extract from the data.
#' @param Xvars a vector of covariate names that will be used in the analysis.
#' 
#' 
#' @details The most common format for storing compositional data is matrix format. The 
#' \code{read.dm} function can then be used to read in the data. Occasionally data
#' will be stored in predator-prey format. In this instance \code{read.pp} can 
#' be used to read in the data. In either function, either an R object or file name 
#' is required to read in the data.
#'
#' @return a matrix where each row represents data for a separate predator-prey combination.
#' 
#' @references Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
#' and Stable Isotope Data, Journal of Statistical Software, In Prep.
#' 
#' @keywords read
#' @examples 
#' # Reading in a diet matrix
#' data(yftDMraw)
#' write.csv(yftDMraw, file = "yftDMraw.csv", row.names = FALSE)
#' yftpp1 <- read.dm(filenm = "yftDMraw.csv", 
#'                   labels = list(FullnessL = "Fullness", DateL = "Date"),
#'                                     Datef = "\%m/\%d/\%Y", diet.ind.start = 12, p = 0.01)
#' data(PreyTaxonSort)
#' val <- apc(x = yftpp1, preyfile = PreyTaxonSort, check = TRUE)
#' node.colsY <- val$cols
#' dietPP <- val$x   # updated diet matrix with "Group" assigned prey taxa codes
#' head(dietPP)
#' 
#' # Reading in a predator-prey matrix
#' data(yftPPraw)
#' write.csv(yftPPraw, file = "yftPPraw.csv", row.names = FALSE)
#' yftpp2 <- read.pp(filenm = "yftPPraw.csv",
#'                   labels = list(PredatorL = "TripSetPredNo", TripSetL = "TripSetNo",
#'                                                   SpeciesL = "Family", FullnessL = "Fullness", 
#'                                                   DateL = "Date", WeightL = "PropW", PreyGrpL = "Family"), 
#'                                                   Datef = "\%m/\%d/\%Y", p = 0.01,
#'                                                   Xvars = c("Lat", "Lon", "Year", "Quarter", "Length", "SST"))
#' data(PreyTaxonSort)
#' pal <- c(topo.colors(10)[1:2], heat.colors(10)[1:2], terrain.colors(25)[1:8])
#' val <- apc(x = yftpp2, preyfile = PreyTaxonSort, palette = pal, check = TRUE)
#' node.colsY <- val$cols
#' dietPP <- val$x   # updated diet matrix with prey taxa codes
#' head(dietPP)
#' @export


#' @rdname read
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



#' @rdname read
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
