plot.diet <- function(x, Xvar, LonID, LatID, mapxlim, mapylim, PredSpID = NULL, SmXvar = NULL,
                      leg.loc = "topleft", col = "black", Factor,
                      database = 'world', PredIDno = NULL, prey.cols = NULL, 
                      filenm = NULL, ...){
 
  if (!inherits(x, "diet"))
    stop("Not a diet object. \n")
  
  if(missing(LonID) | missing(LatID)){
    geoF <- FALSE
    warning("No maps will be produced as co-ordinates have not been supplied.\n")
  }
  else
    geoF <- TRUE
  
  if(missing(Xvar))
    warning("No exploratory analyses of x-variables will be produced as you have not
            specified any x-variables in the Xvar argument.\n")


  if(!(missing(LonID) | missing(LatID))){
    if(missing(mapxlim))
      mapxlim <- range(x[,LonID])
    if(missing(mapylim))
      mapylim <- range(x[,LatID])
  }


  if(!is.null(filenm)){
    # print to screen
     pdf(filenm)
    
  }
    
  expl1 <- expl2 <- expl3 <- expl4 <- expl5 <- expl6 <- expl7 <- expl8 <- NULL
  expl9 <- expl10 <- NULL
 
  if(geoF){
    
    # Mapping
    # map of all the points
    expl1 <- mappoints.data(x[,c(LonID,LatID)], xlim = mapxlim, ylim = mapylim,  
                   database = database)

    # map of all the points with colours according to predator
    if(!is.null(PredSpID)){
      
      expl2 <- MapPredator(x, LonID, LatID, mapxlim, mapylim, database, PredSpID)

      # maps for each prey
      expl3 <- MapPrey(x, LonID, LatID, mapxlim, mapylim, database)

      # maps for each prey by predator
      if(!is.null(PredSpID)){
        expl4 <- MapPreyPredator(x, PredSpID, LonID, LatID, mapxlim, mapylim,
                                 database)
      }
    }
  }
    
     
  
 
  # Summary plots for Xvar
  if(!missing(Xvar)){ 
    # Pairwise plots
     if(length(unique(x[,PredSpID])) > 1)
        Xdat <- x[,c(Xvar,PredSpID)]
     else
        Xdat <- x[,Xvar]
      X.id <- apply(Xdat, 2, function(x) all(is.na(x)))
      Xdat <- Xdat[,!X.id]
      if(any(X.id))
        cat("Columns: ", names(X.id)[X.id], " contains missing values.\n")

      expl5 <- ggpairs(Xdat)

 
      ff <- data.frame(matrix(NA, nrow = 1, ncol = ncol(Xdat)))
      names(ff) <- names(Xdat)
      for(i in 1:ncol(Xdat))
          ff[1,i] <- is.factor(Xdat[,i])
      XvarC <- data.frame(Xdat[,names(ff)[!ff]])
      names(XvarC) <- names(ff)[!ff]
      XvarF <- data.frame(Xdat[,names(ff)[as.matrix(ff)]])
      names(XvarF) <- names(ff)[as.matrix(ff)]

      if(!is.null(XvarC)){
        cat("Exploratory plots for continuous variables are being created.\n")
        # continuous
        p1 <- p2  <-  list()
        bw <- apply(XvarC, 2, function(x) diff(range(x)))/10
        for(i in 1:ncol(XvarC)){
            p1[[i]] <- ggplot(data.frame(x = rep(1, nrow(XvarC)), y = XvarC[,i]), aes(x, y)) +
            geom_boxplot(fill = "light blue")+ xlab("") + ylab(names(XvarC[i]))
            p2[[i]] <- ggplot(data.frame(x = XvarC[,i]), aes(x)) + geom_histogram(binwidth = bw[i]) +
                    ylab("") + xlab(names(XvarC[i]))
        }

        expl6 <- do.call(grid.arrange, p1)
        expl7 <- do.call(grid.arrange, p2)
       }
  }
  

  if(ncol(XvarF) > 0){
  
    cat("Exploratory plots for factors are being plotted.\n")

    if(!is.null(PredSpID)){
       for(i in 1:(ncol(XvarF))){
         expl8 <- FactorPlots(x, XvarF, i, PredSpID)
      }
    }
    else{
      for(i in 1:ncol(XvarF)){         
        x.tab <- data.frame(table(XvarF[,i])/sum(table(XvarF[,i])))
        expl9 <- ggplot(data = x.tab, aes(Freq)) + geom_bar(color = "lightblue") +
          xlab(names(XvarF)[i])
      }
    }
  }

  # Species Distribution and Composition
  plotSpComp <- SpeciesCompBar(x = x, prey.cols = prey.cols, Factor = Factor, Species = "Species")
 
  # Correlation plot
  XvarCF <- x[,Xvar]
  for(i in 1:ncol(XvarCF))
      XvarCF[,i] <- as.numeric(as.vector(XvarCF[,i]))
  Xcor <- cor(as.matrix(XvarCF), method = "spearman", use = "pairwise.complete.obs")
  meltedXcor <- melt(Xcor)
  names(meltedXcor)[3] <- "rho"
  expl10 <- ggplot(data = meltedXcor, aes(x = X1, y = X2, fill = rho)) + 
    geom_tile() + xlab("") +
    ylab("") + ggtitle("Spearman Rank Correlation of X-variables")

 

  # Smooth exploratory plots
  res <- NULL
  
  if(!is.null(SmXvar)){
    if(missing(LatID) | missing(LonID))
      warning("Maps of smoothed variables can not be produced as you have not specified
              LatID and LonID.\n")
    else{
      cat("Note: Producing interpolated maps for SmXvar specified. Note, a GAM is fitted to each 
          SmXvar variable using the mgcv package and assuming normality. Users should check 
          the fit of each model to see whether these assumptions hold.\n\n")
      res <- smplot <- list()
      
      SmXdat <- data.matrix(x[,SmXvar]) 
      palette <- spaMM.colors()
      projection <- "+proj=longlat +datum=WGS84"
      
      for(i in 1:ncol(SmXdat)){
       val <- SmPlots(x, i, SmXvar, SmXdat, LonID, LatID, projection, palette, too.far = 0.05)
       res[[i]] <- val$res
       smplot[[i]] <- val$smplot
      }
      names(res) <- SmXvar
      
    }
    
  }

  
  if(!is.null(filenm)){
    dev.off()
    cat(paste("Plots written to file ", filenm, "\n", sep = ""))
    
  } # Summary Statistics
  sum.x <- summary(x)
  
  
  sum.x2 <- data.frame(nobs = nrow(x), npred = ifelse(is.null(PredIDno), NA, length(unique(x[,PredIDno]))),
                       nprey = length(levels(x$Group)))
  

  
    list(SmGAMOutput = lapply(res,summary), dataS1 = sum.x, dataS2 = sum.x2,
         expl1 = expl1, expl3 = expl3, expl4= expl4, expl5 = expl5,
         expl6 = expl6, expl7 = expl7, expl8 = expl8, expl9 = expl9, expl10 = expl10,
         plotSpComp = plotSpComp, smplot = smplot)
}
