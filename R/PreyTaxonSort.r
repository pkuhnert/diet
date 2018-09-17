#'
#' Prey taxa found in Yellowfin tuna stomach contents
#' 
#' @description Taxonomic groupings of prey found in YFT stomach contents.  
#' 
#' @usage data(PreyTaxonSort)
#' 
#' @format 
#' A data frame with 12 observations on the following 4 variables.
#' \describe{
#' \item{\code{Group}}{Prey groups, identified to Family level of taxonomy, found in YFT stomach contents 
#' with the following levels: Argonautidae, Exocoetidae, Galatheidae, Hemiramphidae, Myctophidae, Nomeidae, Ommastrephidae,
#' Ostraciidae, Phosichthyidae, Portunidae, Scombridae, Tetraodontidae}
#' \item{\code{PreyTaxonSort}}{Abbreviations for prey groups found in YFT stomach contents with the 
#' following levels: A-Gal (Arthropods-Galatheidae), A-Por (Arthropods-Portunidae), F-Exo (Fishes-Exocoetidae),
#' F-Hem (Fishes-Hemiramphidae), F-Myc (Fishes-Myctophidae), F-Nom (Fishes-Nomeidae), F-Ost (Fishes-Ostraciidae),
#' F-Pho (Fishes-Phosichthyidae), F-Sco (Fishes-Scombridae), F-Tet (Fishes-Tetraodontidae), M-Arg 
#' (Molluscs-Argonautidae), M-Omm (Molluscs-Ommastrephidae) }
#' \item{\code{PreyTaxBroad}}{Broad taxonomic grouping of prey found in YFT stomach contents with the 
#' following levels: Arthropods, Fishes, Molluscs}
#' \item{\code{Sort}}{Sorting column that sorts the prey in order for plotting}
#' }
#' 
#' @details 
#' This dataset is used to sort prey by broad taxonomic groupings for output in the 
#' classification tree plots. Abbreviations of prey groups make the classification 
#' tree visually appealing.
#' 
#' @references 
#' Olson RJ, Popp BN, Graham BS, Lopez-Ibarra GA, Galvan-Magana F, Lennert-Cody CE, 
#' Bocanegra-Castillo N, Wallsgrove NJ, Gier E, Alatorre-Ramirez V, Balance LT, Fry B (2010) 
#' Food web inferences of stable isotope spatial patterns in copepods and yellowfin tuna in the 
#' pelagic eastern Pacific Ocean. Prog Oceanogr.
#' 
#' @examples 
#' data(PreyTaxonSort)
#' PreyTaxonSort
"PreyTaxonSort"