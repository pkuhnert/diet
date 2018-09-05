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

 
   
   
   flush.console()
   
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


