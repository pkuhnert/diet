#' @title The Diet Package
#' 
#' @name diet
#' 
#' @description 
#' The diet package builds on two primary R packages for analysing compositional diet data. 
#' The first is a diet tree analysis called \code{dpart}, which takes a predator-prey 
#' response and builds a classification tree (Breiman et al. 1984) using the 
#' recursive partitioning algorithm available in the \code{rpart} R package (Therneau and 
#' Atkinson 1997) and outlined in Kuhnert et al.(2012). The second is an isotope 
#' analysis called igam, that fits a Generalized Additive Model (GAM) using the 
#' \code{mgcv} package developed by Simon Wood (2006). Together, these primary functions 
#' comprise the diet package which allows researchers to analyse trophic 
#' ecology data, form predictions, understand the variation in the data and determine 
#' what predictors, if any, can lead to changes in diet characteristics.
#' 
#' @details 
#'   Package: diet
#'   Type: Package
#'   Version:  1.0
#'   Date:  2018-09-19
#'   License:  GPL
#'   LazyLoad:  yes
#'   
#'   
#' @author Petra Kuhnert and Leanne Duffy
#'   
#' @references 
#'   Breiman, L., Friedman, J.H., Olshen, R.A. and Stone, C.J. (1984) Classification 
#'   and Regression Trees. Wadsworth International.
#'   
#'   Kuhnert, P.M., Duffy, L.M., Young, J.W. and Olson, R.J. (2012) Predicting fish diet
#'   composition using a bagged classification tree approach: a case study using yellowfin
#'   tuna (Thunnus albacares), Marine Biology, 159, 87-100.
#'   
#'   Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
#'   and Stable Isotope Data, Journal of Statistical Software, In Prep. 
#'   
#'   Therneau, T.M. and Atkinson, E.J. (1997) Recursive Partitioning using the RPART 
#'   Routines, Mayo Foundation. 
#'   
#'   Wood, S. (2006) Generalized Additive Models: An introduction with R, Chapman and 
#'   Hall, New York.
#'   
#' @seealso  \code{\link{rpart}}; \code{\link{gam}}
#'   
#' @keywords package
NULL