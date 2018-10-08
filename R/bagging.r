#' Bagging
#' 
#' @description Creates bagged tree estimates from diet data
#' 
#' @param formula a \link{formula}, with a response but no interaction terms 
#' as for the \code{\link{rpart}} function
#' @param data an optional data frame in which to interpret the variables named in the formula
#' @param weights case weights
#' @param subset  optional expression saying that only a subset of the rows of the data should be used in the fit.
#' @param na.action The default action deletes all observations for which \code{y} is missing, but keeps
#' those in which one or more predictors are missing.
#' @param model  if logical: keep a copy of the model frame in the result? If the input value for
#' model is a model frame (likely from an earlier call to the rpart function), then
#' this frame is used rather than constructing new data.
#' @param x keep a copy of the x matrix in the result.
#' @param y keep a copy of the dependent variable in the result. If missing and \code{model}
#' is supplied this defaults to \code{FALSE}.
#' @param parms optional parameters for the splitting function. For classification splitting, 
#' the list can contain any of: the vector of prior probabilities (component prior), the loss 
#' matrix (component loss) or the splitting index (component split). The priors must be 
#' positive and sum to 1. The loss matrix must have zeros on the diagonal and positive 
#' off-diagonal elements. The splitting index can be gini or information. The default 
#' priors are proportional to the data counts, the losses default to 1, and the split 
#' defaults to \code{gini}.
#' @param control options that control details of the \code{rpart} algorithm.
#' @param cost a vector of non-negative costs, one for each variable in the model.
#' Defaults to one for all variables. These are scalings to be applied when
#' considering splits, so the improvement on splitting on a variable is divided by
#' its cost in deciding which split to choose.
#' @param nBaggs numeric. Number of bootstrap samples.
#' @param spatial A list with the following elements:
#' fit  = do spatial bootstrapping 
#' sizeofgrid =  size of spatial tile to sample from (default is 5)
#' nsub = number of sub-samples to take (defaults to no subsampling)
#' ID   = ID in which to subsample from (e.g. TripSetPredNo)
#' (only required if sub-sampling is required)
#' @param Plot plotting the spatial grid with samples (default: no plotting (FALSE))
#' @param predID predator ID
#' @param numCores Number of cores to push the bagging on to. Only available under Unix (default: 1) 
#' @param \dots arguments to be passed to or from other methods.
#' 
#' @usage 
#'   bagging(formula, data, weights, subset, na.action = na.dpart,
#'             model = FALSE, x = FALSE, y = TRUE, parms, control, 
#'             cost, nBaggs,
#'                       spatial = list(fit = FALSE, sizeofgrid = 5, 
#'                       nsub = NULL, ID = NULL, LonID = "Longitude", 
#'                                                LatID = "Latitude"), 
#'                                                Plot = FALSE,
#'                                    predID, numCores = 1, ...)
#'                                                          
#' @details   Users will need to determine whether spatial bootstrapping is required. They can
#' use the \code{\link{resid}} function to examine the residuals from the fit of the
#' model to determine whether this is required.
#' 
#' @return   A list with the following elements:
#' \item{baggs}{tree objects for each \code{B} trees produced.}
#' \item{oob}{numeric vector indicating the samples left as out of bag (oob) samples.}
#' \item{pred.oob}{predicted prey composition for each set of out of bag samples.}
#' \item{pred}{all predicted prey compositions for each bootstrap sample.}
#' \item{resid}{data frame of residuals from the fitted tree for each bootstrap sample.}
#' \item{data}{bootstrap sample dataset}
#' 
#' @references   Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
#' and Stable Isotope Data, Journal of Statistical Software, In Prep. 
#' 
#' Kuhnert PM, Kinsey-Henderson A, Bartley R, Herr A (2010)
#' Incorporating uncertainty in gully erosion calculations using 
#' the random forests modelling approach. Environmetrics 21:493-509. doi:10.1002/env.999
#' 
#' Breiman L (1996) Bagging predictors. Mach Learn 24:123-140. doi:10.1023/A:1018054314350  
#' 
#' Breiman L (1998) Arcing classifiers (with discussion). Ann Stat 26:801-824. doi:10.2307/120055
#' 
#' Breiman L (2001) Random forests. Mach Learn 45:5-32. doi:10.1023/A:1010933404324
#' 
#' 
#' @examples 
#' 
#' # Assigning prey colours for default palette
#' #val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = TRUE)
#' #node.colsY <- val$cols
#' #dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes
#' 
#' 
#' # Bagging
#' # Bagging with NO spatial bootstrapping
#' # N.B. Not run as this takes a while
#' #yft.bag <- bagging(Group ~ Lat + Lon + Year + Quarter + SST  + Length,
#' #                      data = dietPP, weights = W, minsplit = 50,
#' #                       cp = 0.001, nBaggs = 500, predID = "TripSetPredNo")
#' #
#' 
#'                      
#' @import "foreach"
#'  
#' @export
bagging <- function(formula, data, weights, subset, na.action = na.dpart,
      model = FALSE, x = FALSE, y = TRUE, parms, control, cost, nBaggs,
      spatial = list(fit = FALSE, sizeofgrid = 5, nsub = NULL, ID = NULL, LonID = "Longitude", 
                     LatID = "Latitude"), Plot = FALSE,
       predID, numCores = 1, ...){

      # This function produces (naive) bagged estimates using the predator-prey matrix
      # Args:
      # formula   = as for dpart
      # data      = as for dpart (predator-prey df)
      # weights, subset, na.action, method, model, x, y, parms,
      # control, cost   = as for dpart
      # nBaggs = number of bootstraps
      # spatial = list with the following elements:
      #   fit  = do spatial bootstrapping
      #   sizeofgrid =  size of spatial tile to sample from (default is 5)
      #   nsub = number of sub-samples to take (defaults to no subsampling)
      #   ID   = ID in which to subsample from (e.g. TripSetPredNo)
      #          (only required if sub-sampling is required)
      # Plot = plotting the spatial grid with samples (default: no plotting (FALSE))

 

   bSample <- function(dat){
     resample <- function(x, ...) x[sample.int(length(x), ...)]
     tt <- unlist(tapply(1:nrow(dat), as.vector(dat$Group), function(x) resample(x, size = length(x), 
                                                                    replace = TRUE)))

     dat[tt, ]
   } 
   

   m <- match.call(expand.dots = FALSE)
   if(!is.na(match("spatial", names(m))))
      m <- m[-match("spatial", names(m))]
   if(!is.na(match("Plot", names(m))))
      m <- m[-match("Plot", names(m))]
   if(!is.na(match("predID", names(m))))
      m <- m[-match("predID", names(m))]
  if(!is.na(match("numCores", names(m))))
    m <- m[-match("numCores", names(m))]
  

   m$model <- m$method <- m$control <- NULL
   m$x <- m$y <- m$parms <- m$... <- NULL
   m$nBaggs <- NULL
   m$cost <- NULL
   m$na.action <- na.action
   
   m[[1L]] <- as.name("model.frame")
   m <- eval(m, parent.frame())
   Y <- model.extract(m, "response")

    if(missing(weights))
      W <- match.call(expand.dots = FALSE)$weights
   else{
      W <- model.extract(m, "weights")
      data$W <- W
   }
   if(missing(cost)){
     Terms <- attr(m, "terms")
     nvar <- length(attr(Terms, "predvars"))-2
     cost <- rep(1, nvar)
   }
   if(missing(control))
      control <- rpart.control(...)
 
   options(warn = -1)
   
   baggs <- oob <- pred.oob <- pred <- pp <- list()
 
  
   Omat <- formOmat(data, ID = predID)
   resids <- data.frame(matrix(NA, nrow = nrow(Omat), ncol = nBaggs))
   names(resids) <- paste("B", 1:nBaggs, sep = "")
   row.names(resids) <- row.names(Omat)
   
  
   
 #  registerDoMC(numCores)     # only available on unix machines


   #for(i in 1:nBaggs){
    
   result <- foreach(i = 1:nBaggs) %dopar%{
   if(i %% 50 == 0)
       cat("Iteration ", i, "\n")
     
 
    
       if(spatial$fit)
          bsamp <- spatialsamp(x = data, LonID = spatial$LonID, LatID = spatial$LatID,
                               sizeofgrid = spatial$sizeofgrid,
                       ID = spatial$ID, nsub = spatial$nsub, Plot = Plot)
       else
          bsamp <- bSample(data)
       sampid <- as.integer(row.names(bsamp))
       oob.m <- match(as.integer(row.names(data)), as.integer(row.names(bsamp)))
       oob[[i]] <- as.integer(row.names(data))[is.na(oob.m)]
       W.b <- W[sampid]
       baggs[[i]] <- dpart(formula, bsamp, weights = W, na.action = na.action,
              model = model, x = x, y = y, control = control,
              cost = cost, ...)
       pred.oob[[i]] <- predict(baggs[[i]], newdata = data[paste(oob[[i]]),], 
          type = "prob", plot = FALSE)
       pred[[i]] <- predict(baggs[[i]], newdata = data, type = "prob", plot = FALSE)
     ###

     pp[[i]] <- predict(baggs[[i]], newdata = data, type = "prob", pred.type = "predator",
                             predatorID = predID, plot = FALSE)
     d <- Distance(O = Omat[,(ncol(Omat)-ncol(pp[[i]])+1):ncol(Omat)], P = pp[[i]], type = "Hellinger")
     resids[,i] <- d$Dist
     
     
   }
   actualresult <- unlist(result)
  

    Ebag <- list(baggs = baggs, oob = oob, pred.oob = pred.oob, pred = pred, predP = pp,
                resid = resids,  data = data)
    class(Ebag) <- c("bag")
    options(warn = 0)
    
    Ebag
    
}


