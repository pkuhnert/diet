#'
#' Yellowfin Tuna isotope dataset
#' 
#' @description 
#' Yellowfin tuna (YFT) were captured by purse-seine fishing vessels in the eastern
#' tropical Pacific (ETP) between 16 August 2003 and 16 November 2005. Samples of YFT white muscle tissue
#' were taken for stable isotope analysis.
#' 
#' 
#' @format 
#' A data frame with 49 observations on the following 10 variables.
#' \describe{
#' \item{\code{IsoSampNo}}{Unique, nonrepeating number for each isotope sample taken from an individual 
#' YFT or several YFT aggregated into a composite sample }
#' \item{\code{Yr}}{Year in which YFT was captured}
#' \item{\code{Month}}{Month YFT was captured}
#' \item{\code{Date}}{Date YFT was captured}
#' \item{\code{Lat}}{Latitude coordinates where YFT was captured}
#' \item{\code{Lon}}{Longitude coordinates where YFT was captured}
#' \item{\code{SST}}{Sea-surface temperature at location where YFT was captured, measured in degrees C}
#' \item{\code{FL_mm}}{Length, measured as fork length, of YFT in mm}
#' \item{\code{d13C}}{YFT stable isotope ratios relative to the isotope standards 
#' are expressed in the following conventional delta (\eqn{\delta}) notation in parts per thousand (\%):
#' \deqn{\delta X = ([R_{sample}/R_{standard}]-1)\times 1000}
#' where X is \eqn{^{15}N} or \eqn{^{13}C}, and \eqn{R_{sample}} and \eqn{R_{standard}} are the 
#' corresponding ratios of heavy to light isotopes (\eqn{^{15}N/^{14}N} or \eqn{^{14}C/^{13}C}) 
#' in the sample and standard, respectively. 
#' \eqn{R_{standard}} for \eqn{\delta^{15}N} was atmospheric \eqn{N_2} and for \eqn{\delta^{13}C} 
#' was the Peedee Belemnite (PDB) limestone formation.}
#' \item{\code{d15N}}{YFT stable isotope ratios relative to the isotope standards are 
#' expressed in the following conventional delta (\eqn{\delta}) notation in parts per thousand (\%): 
#' \deqn{\delta X = ([R_{sample}/R_{standard}]-1)\times 1000}where X is \eqn{^{15}N} or \eqn{^{13}C}, 
#' and \eqn{R_{sample}} and \eqn{R_{standard}} are the corresponding ratios of heavy to light isotopes
#' (\eqn{^{15}N/^{14}N} or \eqn{^{14}C/^{13}C}) 
#' in the sample and standard, respectively. 
#' \eqn{R_{standard}} for \eqn{\delta^{15}N} was atmospheric \eqn{N_2} and for \eqn{\delta^{13}C} 
#' was the Peedee Belemnite (PDB) limestone formation.}
#' }
#' 
#' @details 
#' YFT were sampled onboard the vessels by observers of the Inter-American Tropical Tuna Commission (IATTC, 2004).  
#' Fish were randomly subsampled from purse-seine sets immediately after capture, and date, position, and 
#' sea-surface temperature (SST) were recorded for each set that yielded samples.  Observers measure
#' fork length of the fish, removed samples of white muscle tissue from the dorsal musculature adjacent
#' to the second dorsal fin and stored them at about -20 degrees C until processed further in the lab.  Subsamples of white muscle 
#' from up to six individual YFT per purse-seine set and size class (<90 and >= 90 cm FL) were combined into one composite sample 
#' for stable isotope analysis. Thus, there were 49 composite samples of 225 fish.
#' 
#' @references 
#' Olson RJ, Popp BN, Graham BS, Lopez-Ibarra GA, Galvan-Magana F, Lennert-Cody CE, 
#' Bocanegra-Castillo N, Wallsgrove NJ, Gier E, Alatorre-Ramirez V, Balance LT, Fry B (2010) 
#' Food web inferences of stable isotope spatial patterns in copepods and yellowfin tuna in the 
#' pelagic eastern Pacific Ocean. Prog Oceanogr.
#' 
"yftISO"