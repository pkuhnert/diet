% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.r
\name{read.dm}
\alias{read.dm}
\alias{read,}
\alias{read.pp}
\title{Read in Diet Matrix}
\usage{
read.dm(object = NULL, filenm = NULL, labels = list(DateL = "Date",
  FullnessL = "Fullness"), Datef = "\%d/\%m/\%Y", diet.ind.start,
  p = 0.01, sp2omit = NULL, predsp = NULL)

read.pp(object = NULL, filenm = NULL, labels = list(PredatorL =
  "PredID", TripSetL = "TripSetNo", SpeciesL = "Species", FullnessL =
  "Fullness", DateL = "Date", WeightL = "Weight", PreyGrpL = "PreyGrp"),
  Datef = "\%d/\%m/\%Y", p = 0.01, Xvars, sp2omit = NULL)
}
\arguments{
\item{object}{matrix. A compositional data matrix containing covariate 
information and diet composition for each predator (optional).}

\item{filenm}{file name. A file name specifying the location of a comma 
delimited file that contains the covariate information and diet composition 
for each predator (optional).}

\item{labels}{a list containing the following: 
\code{DateL} Date column name: default label is "Date".
\code{FullnessL} Fullness column name used to identify full stomachs for analysis: 
default label is "Fullness".
\code{PredatorL} Predator column name: default label is "PredID".
\code{TripSetL} Sampling column name: default label is "TripSetNo".
\code{SpeciesL} Species column name: default label is "Species". This may 
have the same label as PreyGrpL.
\code{WeightL} Weight column name: default label is "Weight". This may 
represent weights, frequencies or binary data.
\code{PreyGrpL} Prey group column name: default label is "PreyGrp". This may 
have the same label as SpeciesL.}

\item{Datef}{date format. See \code{\link{strptime}} for format options.}

\item{diet.ind.start}{integer. Column number where the diet composition data begins.}

\item{p}{numeric. Proportion of prey to omit. If the composition is represented as a 
wet weight, then \code{p} represents the proportion of prey contributing to the wet 
weight consumption of a predator. If the composition is represented as a frequency or 
binary outcome, then \code{p} represents the proportion of prey found in a predator. (default: 0.01)}

\item{sp2omit}{a vector of prey names to omit from the analysis (default: NULL)}

\item{predsp}{Predator species names to extract from the data.}

\item{Xvars}{a vector of covariate names that will be used in the analysis.}
}
\value{
a matrix where each row represents data for a separate predator-prey combination.
}
\description{
Read in diet compositional data in one of two formats. Data 
stored as a compositional matrix can be read in using the \code{read.dm} 
function. Data stored in a predator-prey format can be read in using the 
\code{read.pp} function.
}
\details{
The most common format for storing compositional data is matrix format. The 
\code{read.dm} function can then be used to read in the data. Occasionally data
will be stored in predator-prey format. In this instance \code{read.pp} can 
be used to read in the data. In either function, either an R object or file name 
is required to read in the data.
}
\examples{
# Reading in a diet matrix
#data(yftDMraw)
#write.csv(yftDMraw, file = "yftDMraw.csv", row.names = FALSE)
#yftpp1 <- read.dm(filenm = "yftDMraw.csv", 
#                   labels = list(FullnessL = "Fullness", DateL = "Date"),
#                                     Datef = "\\\%m/\\\%d/\\\%Y", diet.ind.start = 12, p = 0.01)
#data(PreyTaxonSort)
#val <- apc(x = yftpp1, preyfile = PreyTaxonSort, check = TRUE)
#node.colsY <- val$cols
#dietPP <- val$x   # updated diet matrix with "Group" assigned prey taxa codes
#head(dietPP)

# Reading in a predator-prey matrix
#data(yftPPraw)
#write.csv(yftPPraw, file = "yftPPraw.csv", row.names = FALSE)
#yftpp2 <- read.pp(filenm = "yftPPraw.csv",
#                   labels = list(PredatorL = "TripSetPredNo", TripSetL = "TripSetNo",
#                               SpeciesL = "Family", FullnessL = "Fullness", 
#                                DateL = "Date", WeightL = "PropW", PreyGrpL = "Family"), 
#                                 Datef = "\\\%m/\\\%d/\\\%Y", p = 0.01,
#                               Xvars = c("Lat", "Lon", "Year", "Quarter", "Length", "SST"))
#data(PreyTaxonSort)
#pal <- c(topo.colors(10)[1:2], heat.colors(10)[1:2], terrain.colors(25)[1:8])
#val <- apc(x = yftpp2, preyfile = PreyTaxonSort, palette = pal, check = TRUE)
#node.colsY <- val$cols
#dietPP <- val$x   # updated diet matrix with prey taxa codes
#head(dietPP)
}
\references{
Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
and Stable Isotope Data, Journal of Statistical Software, In Prep.
}
\keyword{read}
