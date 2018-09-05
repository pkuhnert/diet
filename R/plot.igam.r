plot.igam <- function(x, Xvar, var.cond = NULL, plot.map = FALSE,
                       database = 'world', leg.pos = "topleft", spgrp = NULL,
                       too.far = 0.1, se.plot = TRUE, zlim = NULL, ...){

  smterms <- unlist(lapply(x$smooth, function(x) x$term))
  m.smterms <- match(Xvar, smterms)
 
  

  if(length(Xvar) == 1 & !any(is.na(m.smterms))){
    # Contribution of a single variable
    m <- unlist(lapply(x$smooth, function(x, Xvar) 
                         match(Xvar, x$term), Xvar))
    select <- (1:length(m))[m & !is.na(m)]
    plot.gam(x, select = select, rug = FALSE, se = TRUE, shade = TRUE)   
    rug(x$data[,Xvar])
    }
  else if(length(Xvar) == 1){
    # Contribution of a single variable which is a factor
    term <- predict(x, type = "terms", se.fit = TRUE)
    term.fit <- data.frame(term$fit)
    term.sefit <- data.frame(term$se.fit)
    # get rid of smooth terms
    smterms <- ncol(term.fit)-length(x$smooth) + 1
    term.fit <- term.fit[,-(smterms:ncol(term.fit))]
    term.sefit <- term.sefit[,-(smterms:ncol(term.sefit))]
    # Get design matrix
    len <- length(x$assign)
    X <- predict.gam(x, type = "lpmatrix")[,1:len]
    if(is.null(x$na.action))
       dat.Xvar <- x$data[,Xvar]
    else
      dat.Xvar <- x$data[,Xvar][-x$na.action]
    predmat <- data.frame(dat.Xvar, terms = X[, 1:len, drop = FALSE] %*% coef(x)[1:len],
                          se.terms = sqrt(rowSums((X[, 1:len, drop = FALSE] %*% x$Vp[1:len, 1:len]) * 
                            X[, 1:len, drop = FALSE])))    
#    predmat <- data.frame(x$data[,Xvar], terms = term.fit[,Xvar], se.terms = term.sefit[,Xvar])
    names(predmat)[1] <- Xvar
    predmat[,Xvar] <- factor(as.vector(predmat[,Xvar]))
    sp.fit <- tapply(predmat$terms, predmat[,Xvar], mean)
    sp.sefit <- tapply(predmat$se.terms, predmat[,Xvar], mean)
    sp <- data.frame(fit = sp.fit, se.fit = sp.sefit)
    sp$lci <- sp$fit - 1.96 * sp$se.fit
    sp$uci <- sp$fit + 1.96 * sp$se.fit
    
    if(!is.null(spgrp)){ 
       sp$group <- spgrp$Group[match(row.names(sp), spgrp[,Xvar])]
       sp$group <- as.factor(sp$group)
       sp$cols <- sp$group
       levels(sp$cols) <- rainbow(length(levels(sp$cols)))
       sp <- sp[order(sp$cols),]
       col.lev <- levels(sp$cols)
       sp$cols <- as.vector(sp$cols)
    }
    else
      sp$cols <- 1
    # Plotting
    par(mar = c(7, 5, 6, 2)+0.1)
    plot(1:nrow(sp), sp$fit, ylim = range(sp$lci, sp$uci), type = "n", 
            axes = FALSE, xlab = "", ylab = paste(x$formula[[2]]))
    arrows(1:nrow(sp), sp$lci, 1:nrow(sp), sp$uci, length = 0.1, angle = 90, code = 3, 
           col = sp$cols)
    points(1:nrow(sp), sp$fit, pch = 15, col = sp$cols)
    axis(side = 2)
    box()
    mtext(side = 1, at = 1:nrow(sp), text = row.names(sp), las = 2, adj = 1, line = 0.5,
          col = as.vector(sp$cols))
    nobs <- table(x$data[,Xvar])
    nlab <- nobs[match(row.names(sp), names(nobs))]  
    mtext(side = 1, at = 1:nrow(sp), text = paste("(n=",nlab,")", sep = ""), adj = 1, line = 3.5,
          cex = 0.8, las = 2)
  
    if(!is.null(spgrp))
      legend("topright", col = col.lev, lty = 1, legend = levels(sp$group), 
           bty = "n", cex = 0.8)
    title(main = Xvar)
    par(mar = c(5,5,6,2)+0.1)
    
  }
  else if(length(Xvar) == 2 & is.null(var.cond)){
   #2D plot
  
    if(is.null(zlim)){
      zlim <- range(x$y, na.rm = TRUE)
      miny <- min(x$y)
      maxy <- max(x$y)
    }
    else{
      miny <- zlim[1]
      maxy <- zlim[2]
    }
    par(fig = c(0,1,0,1), mar = c(5,4,4,2)+0.1)
    vis(x, Xvar, too.far = too.far, color = "topo", plot.type = "contour", 
            zlim = zlim, contour.col = "slategray", main = "", se.plot = se.plot)
    points(x$data[,Xvar[1]], x$data[,Xvar[2]], pch = 16)
    if(plot.map)
       map(database, add = TRUE, col = "black", fill = FALSE)
    box()
  
    # legend
    par(fig = c(0.2, 0.8, 0.9, 1.0), mar = c(2, 4, 0.5, 2) + 0.1, new = TRUE)
    z <- matrix(seq(0, 1, length = 20), ncol = 1)
    image(1:20, 1, z, col = topo.colors(20),
          xlab = "", ylab = "", axes = FALSE)
    axis(side = 1, round(seq(miny, maxy, length = 20), 2), at = 1:20, cex.axis = 0.8)
    box()
  
    par(fig = c(0,1,0,1), mar = c(5,4, 4, 2)+0.1)
  }
  else if(length(Xvar) == 2 & !is.null(var.cond)){
    #2D plot
    if(is.null(zlim)){
      zlim <- range(x$y, na.rm = TRUE)
      miny <- min(x$y)
      maxy <- max(x$y)
    }
    else{
      miny <- zlim[1]
      maxy <- zlim[2]
    }
    
    # Conditioning on a variable
    if(!is.factor(x$data[,var.cond]))
      stop(paste("Cannot produce a plot conditioned by ", var.cond, " as ", var.cond, 
           " is not a factor. Please refit the model with ",var.cond, " as a factor.\n", sep = ""))
     nm <- levels(x$data[,var.cond])
    if(!is.null(spgrp)){
      group <- spgrp$Group[match(nm, spgrp[,var.cond])]
      nm <- nm[order(group)]
      group <- sort(group)
      group.cols <- group
      levels(group.cols) <- rainbow(length(levels(group)))
      group.cols <- as.vector(group.cols)      
  }
    id <- match(var.cond, names(x$data))
   # do 4 plots per page
    nplots <- ceiling(length(levels(x$data[,var.cond]))/4)
    k <- 0
    for(i in 1:nplots){
      sp <- (4*i - 3):(4*i)
      if(length(levels(x$data[,var.cond])) < sp[4])
         np <- (1:4)[match(length(levels(x$data[,var.cond])), sp)]      
      else
        np <- 4
      for(j in 1:np){
       
        k <- k+1
        # margins
        if(j == 1)
          par(new = FALSE, fig = c(0, 0.5, 0.45, 0.9), mar = c(0, 4, 4, 0))
        else if(j == 2)
          par(new = TRUE, fig = c(0.5, 1.0, 0.45, 0.9), mar = c(0, 0, 4, 4))
        else if(j == 3)
          par(new = TRUE, fig = c(0, 0.5, 0, 0.45), mar = c(4, 4, 0, 0), ask = TRUE)
        else
          par(new = TRUE, fig = c(0.5, 1.0, 0, 0.45), mar = c(4, 0, 0, 4))
        varC <- vector(mode = "list", length = 1)
        names(varC) <- var.cond
        varC[[1]] <- nm[k]
       
        vis(x, view = Xvar, cond = varC, plot.type = "contour", xlab = "", ylab = "",
                color = "topo", type = "response", main = "", too.far = too.far,
                zlim = zlim, contour.col = "slategray", axes = FALSE, se.plot = se.plot)
        
        points(x$data[,Xvar[1]][x$data[,id] == nm[k]], 
               x$data[,Xvar[2]][x$data[,id] == nm[k]], pch = 16)
        if(plot.map)
          map(database, add = TRUE, col = "gold3", fill = TRUE)
        box()
        
        # axes
        if(j == 1)
          axis(side = 3)
        else if(j == 2)
          axis(side = 4)
        else if(j == 3)
          axis(side = 2)
        else
          axis(side = 1)
        if(!is.null(spgrp))
          legend(leg.pos, legend = nm[k], bty = "n", text.col = group.cols[k])
        else
          legend(leg.pos, legend = nm[k], bty = "n")
        # Titles
        
        # legend
        par(fig = c(0.2, 0.8, 0.85, 1), mar = c(3, 4, 1, 2) + 0.1, new = TRUE)
        z <- matrix(seq(0, 1, length = 20), ncol = 1)
        image(1:20, 1, z, col = topo.colors(20),
              xlab = "", ylab = "", axes = FALSE, zlim = zlim)
        axis(side = 1, round(seq(miny, maxy, length = 20), 2), at = 1:20, cex.axis = 0.8)
        box()
        
        par(fig = c(0,1,0,1), mar = c(5,4, 4, 2)+0.1, ask = FALSE)
        
      } 
      if(k <=2){
        par(new = TRUE, fig = c(0, 1.0, 0.35, 0.9), mar = c(5, 4, 4, 2) + 0.1)
        plot(0, 1, xlim = c(0, 1), ylim = c(0,1), axes = FALSE, xlab = "",
             ylab = "", type = "n")
        mtext(side = 2, paste(Xvar[2]), at = 0.5, line = 2.5)
        mtext(side = 1, paste(Xvar[1]), at = 0.5, line = 3)
      }
      else{          
        par(new = TRUE, fig = c(0, 1, 0, 1), mar = c(5, 4, 4, 2) + 0.1)
        plot(0, 1, xlim = c(0, 1), ylim = c(0,1), axes = FALSE, xlab = "",
             ylab = "", type = "n")
        mtext(side = 2, paste(Xvar[2]), at = 0.5, line = 2.5)
        mtext(side = 1, paste(Xvar[1]), at = 0.5, line = 3)
      }
    }
    
    
    

     
      }
      
          
  
  else stop("The combination you are specifying is not an option in this function.")
  
  par(fig = c(0, 1, 0, 1), mar = c(5, 5, 6, 2)+0.1)
  
  
}