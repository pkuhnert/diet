\name{yftDMraw}
\alias{yftDMraw}
\docType{data}
\title{Yellowfin Tuna Stomach Contents Matrix Raw Data 
%%   ~~ data name/kind ... ~~
}
\description{Yellowfin tuna (YFT) were captured by purse-seine fishing vessels in the eastern
tropical Pacific (ETP) between 16 August 2003 and 16 November 2005. 
}
\usage{data(yftDMraw)}
\format{
  A data frame with 571 observations on the following 40 variables.
  \describe{
    \item{\code{TripSetNo}}{A unique number representing the purse-seine vessel's Trip and Set number}
    \item{\code{TripSetPredNo}}{A unique number representing the purse-seine vessel's Trip and Set number, and the YFT sample number}
    \item{\code{Date}}{Date YFT was captured}
    \item{\code{Quarter}}{Quarter in which YFT was captured}
    \item{\code{Year}}{Year in which YFT was captured}
    \item{\code{Lat}}{Latitude coordinates where YFT was captured}
    \item{\code{Lon}}{Longitude coordinates where YFT was captured}
    \item{\code{SST}}{Sea-surface temperature at location where YFT was captured, measured in degrees C}
    \item{\code{Length}}{Length, measured as fork length, of YFT in mm}
    \item{\code{PredSpp}}{Predator Species Name}
    \item{\code{Fullness}}{Visual estimate of YFT stomach fullness, expressed as a percentage}
    \item{\code{Argonautidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Balistidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{BivalviaC}}{A prey item, identified to Class level of taxonomy, found in YFT stomach contents}
    \item{\code{Bramidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Carangidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Clupeidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Coryphaenidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{DecapodaO}}{A prey item, identified to Order level of taxonomy, found in YFT stomach contents}
    \item{\code{Engraulidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Exocoetidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Galatheidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{GastropodaC}}{A prey item, identified to Class level of taxonomy, found in YFT stomach contents}
    \item{\code{Gempylidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Hemiramphidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Hemisquillidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Loliginidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Myctophidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Nomeidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Octopodidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Ommastrephidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Ostraciidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Penaeidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Phosichthyidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Portunidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Scombridae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Squillidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{StomatopodaO}}{A prey item, identified to Order level of taxonomy, found in YFT stomach contents}
    \item{\code{Tetraodontidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
    \item{\code{Thysanoteuthidae}}{A prey item, identified to Family level of taxonomy, found in YFT stomach contents}
  }
}
\details{
Yellowfin Tuna were sampled onboard
the vessels by observers of the Inter-American Tropical Tuna Commission (IATTC, 2004).  Fish
were randomly subsampled from purse-seine sets immediately after capture, and date, position, and 
sea-surface temperature (SST) were recorded for each set that yielded samples.  Observers measured
fork length of the fish, removed samples of white muscle tissue from the dorsal musculature adjacent
to the second dorsal fin and stored them at about -20 degrees C until processed further in the lab. 
Data is in matrix format where each row represents one YFT.

}
\source{Olson RJ, Popp BN, Graham BS, Lopez-Ibarra GA, Galvan-Magana F, Lennert-Cody CE, 
Bocanegra-Castillo N, Wallsgrove NJ, Gier E, Alatorre-Ramirez V, Balance LT, Fry B (2010) 
Food web inferences of stable isotope spatial patterns in copepods and yellowfin tuna in the 
pelagic eastern Pacific Ocean. Prog Oceanogr.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{Olson RJ, Popp BN, Graham BS, Lopez-Ibarra GA, Galvan-Magana F, Lennert-Cody CE, 
Bocanegra-Castillo N, Wallsgrove NJ, Gier E, Alatorre-Ramirez V, Balance LT, Fry B (2010) 
Food web inferences of stable isotope spatial patterns in copepods and yellowfin tuna in the 
pelagic eastern Pacific Ocean. Prog Oceanogr.

IATTC, 2004. Annual Report of the Inter-American Tropical Tuna Commission, 2003. Inter-American
Tropical Tuna Commission, pp. 98.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(yftDMraw)
}
\keyword{datasets}
