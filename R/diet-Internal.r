#' Internal functions used in the diet package
#' 
#' @description These are internal functions used in the \code{diet} packagethat are not intended to 
#' be directly called by the user.
#' 
#' @keywords internal
#' 
#' @import ggplot2
#' 

#' @rdname diet-Internal
rpart.matrix <- function (frame) 
{
  if (!inherits(frame, "data.frame")) 
    return(as.matrix(frame))
  frame$"(weights)" <- NULL
  terms <- attr(frame, "terms")
  if (is.null(terms)) 
    predictors <- names(frame)
  else {
    a <- attributes(terms)
    predictors <- as.character(a$variables)[-1L]
    removals <- NULL
    if ((TT <- a$response) > 0L) {
      removals <- TT
      frame[[predictors[TT]]] <- NULL
    }
    if (!is.null(TT <- a$offset)) {
      removals <- c(removals, TT)
      frame[[predictors[TT]]] <- NULL
    }
    if (!is.null(removals)) 
      predictors <- predictors[-removals]
    labels <- a$term.labels
    if (abs(length(labels) - length(predictors)) > 0) 
      predictors <- predictors[match(labels, predictors)]
  }
  factors <- sapply(frame, function(x) !is.null(levels(x)))
  characters <- sapply(frame, is.character)
  if (any(factors | characters)) {
    for (preds in predictors[characters]) frame[[preds]] <- as.factor(frame[[preds]])
    factors <- factors | characters
    column.levels <- lapply(frame[factors], levels)
    for (preds in predictors[factors]) frame[[preds]] <- as.numeric(frame[[preds]])
    x <- as.matrix(frame)
    attr(x, "column.levels") <- column.levels
  }
  else x <- as.matrix(frame[predictors])
  class(x) <- "rpart.matrix"
  x
}


#' @rdname diet-Internal
#' @importFrom "grDevices" "dev.cur"
rpart.branch <- function (x, y, node, branch)
  {
    
    if (missing(branch)) {
      if (exists(parms <- paste(".rpart.parms", dev.cur(),
                                sep = "."), envir = .GlobalEnv)) {
        parms <- get(parms, envir = .GlobalEnv)
        branch <- parms$branch
      }
      else branch <- 0
    }
    
    is.left <- (node%%2 == 0)
    node.left <- node[is.left]
    parent <- match(node.left/2, node)
    sibling <- match(node.left + 1, node)
    temp <- (x[sibling] - x[is.left]) * (1 - branch)/2
    xx <- rbind(x[is.left], x[is.left] + temp, x[sibling] - temp,
                x[sibling], NA)
    yy <- rbind(y[is.left], y[parent], y[parent], y[sibling],
                NA)
    list(x = xx, y = yy, nodeL = node[is.left], nodeR = node[sibling])
  }



#' @rdname diet-Internal
rpconvert <- function (x) 
{
  if (!inherits(x, "rpart")) 
    stop("x does not appear to be an rpart object")
  ff <- x$frame
  if (is.null(ff$splits)) {
    warning("x not converted")
    return(x)
  }
  ff$splits <- NULL
  ff$wt <- ff$n
  xlev <- attr(x, "xlevels")
  if (length(xlev) > 0) {
    zz <- as.numeric(names(xlev))
    names(xlev) <- attr(x$terms, "term.labels")[zz]
    attr(x, "xlevels") <- xlev
  }
  if (x$method == "class") {
    temp <- cbind(ff$yval, ff$yval2, ff$yprob)
    dimnames(temp) <- NULL
    ff$yval2 <- temp
    ff$yprob <- NULL
    x$frame <- ff
    temp <- rpart:::rpart.class(c(1, 1, 2, 2), NULL, wt = c(1, 1, 
                                                            1, 1))
    #temp <- rpart.class(c(1, 1, 2, 2), NULL, wt = c(1, 1, 
    #                                                        1, 1))
    x$functions <- list(summary = temp$summary, print = temp$print, 
                        text = temp$text)
  }
  else if (x$method == "anova") {
    x$frame <- ff
    temp <- rpart:::rpart.anova(1:5, NULL, wt = rep(1, 5))
   # temp <- rpart.anova(1:5, NULL, wt = rep(1, 5))
    x$functions <- list(summary = temp$summary, text = temp$text)
  }
  else {
    ff$yval2 <- cbind(ff$yval, ff$yval2)
    x$frame <- ff
    temp <- rpart:::rpart.poisson(1:5, NULL, wt = rep(1, 5))
   # temp <- rpart.poisson(1:5, NULL, wt = rep(1, 5))
    x$functions <- list(summary = temp$summary, text = temp$text)
  }
  class(x) <- "rpart"
  x
}




#' @rdname diet-Internal
#' @import "lattice"
SpeciesCompBar <- function(x, prey.cols, Factor, Species){
  
  ######################################################  
  # Overall Summary of data  
  ########################################################
  
  # Distribution
  prey.tab <- table(x$Group)/sum(table(x$Group))
  
  bcDist <- barchart(prey.tab ~ names(prey.tab),  scales = list(x = list(rot = 45)), 
                     main = paste("Distribution of ", Species), col = "skyblue", ylab = "Proportion")
  #plot(val)
  
  # Composition
  unF <- levels(x[,Factor])
  x.tab <- list()
  for(j in 1:length(unF)){
    datF <- x[x[,Factor] == unF[j],]
    x.tab[[j]] <- tapply(datF$W, datF$Group, sum)
    x.tab[[j]][is.na(x.tab[[j]])] <- 0
  }
  names(x.tab) <- unF
  
  
  n <- unlist(lapply(x.tab, length))
  
  
  x.tab.df <- data.frame(cbind(unlist(x.tab), rep(as.vector(levels(x$Group)), length(levels(x[,Factor]))), 
                               rep(names(n), n)))
  
  names(x.tab.df) <- c("y", "Factor", "Group")
  row.names(x.tab.df) <- NULL
  x.tab.df$y <- as.numeric(as.vector(x.tab.df$y))
  
  tab <- tapply(x.tab.df$y, x.tab.df$Group, function(x) x/sum(x))
  x.bp <- data.frame(cbind(unlist(tab), rep(names(tab), lapply(tab, length)),
                           rep(levels(x$Group), length(tab))))
  names(x.bp) <- c("y", "Factor", "Group")
  x.bp$y <- as.numeric(as.vector(x.bp$y))
  x.bp$y[is.nan(x.bp$y)] <- 0
  if(!is.null(prey.cols))
    x.bp$Group <- factor(as.vector(x.bp$Group), levels = names(prey.cols))
  
  
  if(!is.null(prey.cols)){
    trel.def <- trellis.par.get("superpose.polygon")
    trel.def$col <- prey.cols
    trellis.par.set("superpose.polygon", trel.def)
  }
  
  bc <- barchart(y ~ Factor, groups = x.bp$Group,  data = x.bp, stack = TRUE, 
                 scales = list(x = list(rot = 45)),
                 auto.key = list(space = "right"), ylim = c(0,1), 
                 main = paste("Composition of", Species), ylab = "Proportion")
  #plot(val)
  
  ######################################################  
  # Summary of data by year 
  ########################################################
  
  unyr <- sort(unique(x$Year))
  bc.yr <- list()
  for(i in 1:length(unyr)){
    dat <- x[x$Year == unyr[i] & !is.na(x$Year),]
    prey.tab <- table(dat$Group)/sum(table(dat$Group))
    bc.yr[[i]] <- barchart(prey.tab ~ names(prey.tab),  scales = list(x = list(rot = 45)), 
                           main = paste("Distribution of ", Species, ": ", unyr[i]), col = "skyblue", ylim = c(0,1), ylab = "Proportion")
    
  }
  
  for(i in 1:length(unyr)){
    dat <- x[x$Year == unyr[i] & !is.na(x$Year),]
    
    unF <- levels(dat[,Factor])
    x.tab <- list()
    for(j in 1:length(unF)){
      datF <- dat[dat[,Factor] == unF[j],]
      x.tab[[j]] <- tapply(datF$W, datF$Group, sum)
      x.tab[[j]][is.na(x.tab[[j]])] <- 0
    }
    names(x.tab) <- unF
    
    
    n <- unlist(lapply(x.tab, length))
    
    
    x.tab.df <- data.frame(cbind(unlist(x.tab), rep(as.vector(levels(dat$Group)), length(levels(dat[,Factor]))), 
                                 rep(names(n), n)))
    
    names(x.tab.df) <- c("y", "Factor", "Group")
    row.names(x.tab.df) <- NULL
    x.tab.df$y <- as.numeric(as.vector(x.tab.df$y))
    
    tab <- tapply(x.tab.df$y, x.tab.df$Group, function(x) x/sum(x))
    x.bp <- data.frame(cbind(unlist(tab), rep(names(tab), lapply(tab, length)),
                             rep(levels(dat$Group), length(tab))))
    names(x.bp) <- c("y", "Factor", "Group")
    x.bp$y <- as.numeric(as.vector(x.bp$y))
    x.bp$y[is.nan(x.bp$y)] <- 0
    if(!is.null(prey.cols))
      x.bp$Group <- factor(as.vector(x.bp$Group), levels = names(prey.cols))
    
    
    if(!is.null(prey.cols)){
      trel.def <- trellis.par.get("superpose.polygon")
      trel.def$col <- prey.cols
      trellis.par.set("superpose.polygon", trel.def)
    }
    
    bc.grp <- barchart(y ~ Factor, groups = x.bp$Group,  data = x.bp, stack = TRUE, 
                       scales = list(x = list(rot = 45)),
                       auto.key = list(space = "right"), ylim = c(0,1), 
                       main = paste("Composition of", Species, ": ", unyr[i]), ylab = "Proportion")
    
  }
  
  list(bcDist = bcDist, bc = bc, bc.yr = bc.yr, bc.grp = bc.grp)
  
}


#' @rdname diet-Internal
subsample <- function(dat, ID, n){

x <- NULL
unID <- unique(ID)
t.ID <- table(ID)
if(all(t.ID <= n))
  x <- dat
else{
  # do subsampling
  for(i in 1:length(unID)){
    subdat <- subset(dat, ID == unID[i])
    if(is.null(n))
      x <- rbind(x, subdat[sample(1:nrow(subdat), size = nrow(subdat), replace = TRUE),])
    else
      if(nrow(subdat) <= n)
        x <- rbind(x, subdat[sample(1:nrow(subdat), size = nrow(subdat), replace = TRUE),])
      else
        x <- rbind(x, subdat[sample(1:nrow(subdat), size = n, replace = TRUE),])
  }
}



x

}




#' @rdname diet-Internal
#' @import "graphics"
spatialsamp <- function(x, LonID, LatID, sizeofgrid = 5, ID = NULL, nsub = NULL, Plot = FALSE){


lon.seq <- seq(min(x[,LonID], na.rm = TRUE), max(x[,LonID], na.rm = TRUE),
               by = sizeofgrid)
lat.seq <- seq(from = min(x[,LatID], na.rm = TRUE), to = max(x[,LatID], na.rm = TRUE),
               by = sizeofgrid)

if(Plot){
  plot(x[,LonID], x[,LatID], pch = 16, cex = 0.7)
  abline(h = lat.seq, lty = 3, col = "grey")
  abline(v = lon.seq, lty = 3, col = "grey")
}

gx <- cut(x[,LonID], lon.seq, include.lowest = TRUE)
gy <- cut(x[,LatID], lat.seq, include.lowest = TRUE)
ugx <- unique(gx)
ugy <- unique(gy)
subsamp <- NULL
k <- 1

for(i in 1:length(ugx)){
  for(j in 1:length(ugy)){
    
    ids <- (1:nrow(x))[gx %in% ugx[i] & gy %in% ugy[j]]
    samp <- x[ids,]
    if(Plot)
      with(samp, points(Lon, Lat, pch = 16, col = k, cex = 0.7))
    
    if(is.null(nsub)){
      Ssamp <- sample(1:nrow(samp), nrow(samp), replace = TRUE)
      subsamp <- rbind(subsamp, samp[Ssamp,])
    }
    else{
      Ssamp <- subsample(dat = samp, ID = samp[,ID], n = nsub)
      subsamp <- rbind(subsamp, Ssamp)
    }
    k <- k+1
  }
}

subsamp

}






#' @rdname diet-Internal
formOmat <- function(object, ID){
  
  
  id <- object[,ID]
  
  Omat <- matrix(0, nrow = length(unique(id)), ncol = length(levels(object$Group))+1)
  Omat <- data.frame(Omat)
  names(Omat) <- c(ID, levels(object$Group))
  Omat[,ID] <- unique(id)
  
  
  Xvars <- NULL
  for(i in 1:length(unique(id))){
    sub <- object[id == id[i],]
    Xvars <- rbind(Xvars, sub[,-c(ncol(sub), ncol(sub)-1)])
    Omat[Omat[,ID] == id[i],][,as.vector(sub$Group)] <- sub$W
  }
  
  id2 <- match(Omat[,ID], id)
  tmp <- object[id2,][,-c(match(ID, names(object)), ncol(object), ncol(object)-1)]
  
  OmatN <- cbind(Omat[,1], tmp, Omat[,-1])
  names(OmatN) <- c(names(Omat)[1], names(tmp), names(Omat)[-1])
  OmatN
}

#' @rdname diet-Internal
#' @import "grDevices"
#' @importFrom "mgcv" "gam"
#' @importFrom "mgcv" "predict.gam"
#' @importFrom "mgcv" "exclude.too.far"
Vis.Gam <- function (x, view = NULL, cond = list(), n.grid = 30, too.far = 0, 
                     col = NA, color = "heat", contour.col = NULL, se = -1, type = "link", 
                     plot.type = "persp", zlim = NULL, nCol = 50, ...) 
{
  fac.seq <- function(fac, n.grid) {
    fn <- length(levels(fac))
    gn <- n.grid
    if (fn > gn) 
      mf <- factor(levels(fac))[1:gn]
    else {
      ln <- floor(gn/fn)
      mf <- rep(levels(fac)[fn], gn)
      mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
      mf <- factor(mf, levels = levels(fac))
    }
    mf
  }
  dnm <- names(list(...))
  v.names <- names(x$var.summary)
  if (is.null(view)) {
    k <- 0
    view <- rep("", 2)
    for (i in 1:length(v.names)) {
      ok <- TRUE
      if (is.matrix(x$var.summary[[i]])) 
        ok <- FALSE
      else if (is.factor(x$var.summary[[i]])) {
        if (length(levels(x$var.summary[[i]])) <= 1) 
          ok <- FALSE
      }
      else {
        if (length(unique(x$var.summary[[i]])) == 1) 
          ok <- FALSE
      }
      if (ok) {
        k <- k + 1
        view[k] <- v.names[i]
      }
      if (k == 2) 
        break
    }
    if (k < 2) 
      stop("Model does not seem to have enough terms to do anything useful")
  }
  else {
    if (sum(view %in% v.names) != 2) 
      stop(gettextf("view variables must be one of %s", 
                    paste(v.names, collapse = ", ")))
    for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], 
                                 c("numeric", "factor"))) 
      stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
  }
  ok <- TRUE
  for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
    if (length(levels(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  else {
    if (length(unique(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  if (!ok) 
    stop(gettextf("View variables must contain more than one value. view = c(%s,%s).", 
                  view[1], view[2]))
  if (is.factor(x$var.summary[[view[1]]])) 
    m1 <- fac.seq(x$var.summary[[view[1]]], n.grid)
  else {
    r1 <- range(x$var.summary[[view[1]]])
    m1 <- seq(r1[1], r1[2], length = n.grid)
  }
  if (is.factor(x$var.summary[[view[2]]])) 
    m2 <- fac.seq(x$var.summary[[view[2]]], n.grid)
  else {
    r2 <- range(x$var.summary[[view[2]]])
    m2 <- seq(r2[1], r2[2], length = n.grid)
  }
  v1 <- rep(m1, n.grid)
  v2 <- rep(m2, rep(n.grid, n.grid))
  newd <- data.frame(matrix(0, n.grid * n.grid, 0))
  for (i in 1:length(x$var.summary)) {
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) {
      ma <- x$var.summary[[i]]
      if (is.numeric(ma)) 
        ma <- ma[2]
    }
    if (is.matrix(x$var.summary[[i]])) 
      newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(x$var.summary[[i]]), 
                          byrow = TRUE)
    else newd[[i]] <- rep(ma, n.grid * n.grid)
  }
  names(newd) <- v.names
  newd[[view[1]]] <- v1
  newd[[view[2]]] <- v2
  if (type == "link") 
    zlab <- paste("linear predictor")
  else if (type == "response") 
    zlab <- type
  else stop("type must be \"link\" or \"response\"")
  fv <- predict.gam(x, newdata = newd, se.fit = TRUE, type = type)
  z <- fv$fit
  if (too.far > 0) {
    ex.tf <- exclude.too.far(v1, v2, x$model[, view[1]], 
                             x$model[, view[2]], dist = too.far)
    fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
  }
  if (is.factor(m1)) {
    m1 <- as.numeric(m1)
    m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
  }
  if (is.factor(m2)) {
    m2 <- as.numeric(m2)
    m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
  }
  if (se <= 0) {
    old.warn <- options(warn = -1)
    av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), n.grid, 
                 n.grid - 1)
    options(old.warn)
    max.z <- max(z, na.rm = TRUE)
    z[is.na(z)] <- max.z * 10000
    z <- matrix(z, n.grid, n.grid)
    surf.col <- t(av) %*% z %*% av
    surf.col[surf.col > max.z * 2] <- NA
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      min.z <- min(fv$fit, na.rm = TRUE)
      max.z <- max(fv$fit, na.rm = TRUE)
    }
    surf.col <- surf.col - min.z
    surf.col <- surf.col/(max.z - min.z)
    surf.col <- round(surf.col * nCol)
    con.col <- 1
    if (color == "heat") {
      pal <- heat.colors(nCol)
      con.col <- 3
    }
    else if (color == "topo") {
      pal <- topo.colors(nCol)
      con.col <- 2
    }
    else if (color == "cm") {
      pal <- cm.colors(nCol)
      con.col <- 1
    }
    else if (color == "terrain") {
      pal <- terrain.colors(nCol)
      con.col <- 2
    }
    else if (color == "gray" || color == "bw") {
      pal <- gray(seq(0.1, 0.9, length = nCol))
      con.col <- 1
    }
    else stop("color scheme not recognised")
    if (is.null(contour.col)) 
      contour.col <- con.col
    surf.col[surf.col < 1] <- 1
    surf.col[surf.col > nCol] <- nCol
    if (is.na(col)) 
      col <- pal[as.array(surf.col)]
    z <- matrix(fv$fit, n.grid, n.grid)
    if (plot.type == "contour") {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",main=zlab"), ",...)", 
                    sep = "")
      if (color != "bw") {
        txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)", 
                     ifelse("add" %in% dnm, "", ",add=TRUE"), ",...)", 
                     sep = "")
      }
      else {
        txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)", 
                     stub, sep = "")
      }
    }
    else {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("zlab" %in% dnm, "", ",zlab=zlab"), ",...)", 
                    sep = "")
      if (color == "bw") {
        txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", 
                     stub, sep = "")
      }
      else {
        txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                     stub, sep = "")
      }
    }
  }
  else {
    if (color == "bw" || color == "gray") {
      subs <- paste("grey are +/-", se, "s.e.")
      lo.col <- "gray"
      hi.col <- "gray"
    }
    else {
      subs <- paste("red/green are +/-", se, "s.e.")
      lo.col <- "green"
      hi.col <- "red"
    }
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      max.z <- max(fv$fit + fv$se.fit * se, na.rm = TRUE)
      min.z <- min(fv$fit - fv$se.fit * se, na.rm = TRUE)
      zlim <- c(min.z, max.z)
    }
    z <- fv$fit - fv$se.fit * se
    z <- matrix(z, n.grid, n.grid)
    if (plot.type == "contour") 
      warning("sorry no option for contouring with errors: try plot.gam")
    stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                  ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
                                                                         dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, 
                                                                                                        "", ",sub=subs"), ",...)", sep = "")
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=lo.col"), stub, sep = "")
    z <- fv$fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=\"black\""), stub, sep = "")
    z <- fv$fit + se * fv$se.fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=hi.col"), stub, sep = "")
  }
  mat <- expand.grid(x = m1, y = m2)
  mat$z <- as.vector(z)
  
  list(mat = mat, x = m1, y = m2, z = z)
}


#' @rdname diet-Internal
vis <- function(x, view = NULL, cond = list(), n.grid = 30, too.far = 0, col = NA, color = "heat", contour.col = NULL, se = -1, type = "link", 
                plot.type = "persp", zlim = NULL, nCol = 50, se.plot = TRUE, ...) UseMethod("vis")

#' @rdname diet-Internal
vis.igam <- function (x, view = NULL, cond = list(), n.grid = 30, too.far = 0, 
                      col = NA, color = "heat", contour.col = NULL, se = -1, type = "link", 
                      plot.type = "persp", zlim = NULL, nCol = 50, se.plot = TRUE, ...) 
{
  
  fac.seq <- function(fac, n.grid) {
    fn <- length(levels(fac))
    gn <- n.grid
    if (fn > gn) 
      mf <- factor(levels(fac))[1:gn]
    else {
      ln <- floor(gn/fn)
      mf <- rep(levels(fac)[fn], gn)
      mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
      mf <- factor(mf, levels = levels(fac))
    }
    mf
  }
  dnm <- names(list(...))
  v.names <- names(x$var.summary)
  if (is.null(view)) {
    k <- 0
    view <- rep("", 2)
    for (i in 1:length(v.names)) {
      ok <- TRUE
      if (is.matrix(x$var.summary[[i]])) 
        ok <- FALSE
      else if (is.factor(x$var.summary[[i]])) {
        if (length(levels(x$var.summary[[i]])) <= 1) 
          ok <- FALSE
      }
      else {
        if (length(unique(x$var.summary[[i]])) == 1) 
          ok <- FALSE
      }
      if (ok) {
        k <- k + 1
        view[k] <- v.names[i]
      }
      if (k == 2) 
        break
    }
    if (k < 2) 
      stop("Model does not seem to have enough terms to do anything useful")
  }
  else {
    if (sum(view %in% v.names) != 2) 
      stop(paste(c("view variables must be one of", v.names), 
                 collapse = ", "))
    for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], 
                                 c("numeric", "factor"))) 
      stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
  }
  ok <- TRUE
  for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
    if (length(levels(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  else {
    if (length(unique(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  if (!ok) 
    stop(paste("View variables must contain more than one value. view = c(", 
               view[1], ",", view[2], ").", sep = ""))
  if (is.factor(x$var.summary[[view[1]]])) 
    m1 <- fac.seq(x$var.summary[[view[1]]], n.grid)
  else {
    r1 <- range(x$var.summary[[view[1]]])
    m1 <- seq(r1[1], r1[2], length = n.grid)
  }
  if (is.factor(x$var.summary[[view[2]]])) 
    m2 <- fac.seq(x$var.summary[[view[2]]], n.grid)
  else {
    r2 <- range(x$var.summary[[view[2]]])
    m2 <- seq(r2[1], r2[2], length = n.grid)
  }
  v1 <- rep(m1, n.grid)
  v2 <- rep(m2, rep(n.grid, n.grid))
  newd <- data.frame(matrix(0, n.grid * n.grid, 0))
  for (i in 1:length(x$var.summary)) {
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) {
      ma <- x$var.summary[[i]]
      if (is.numeric(ma)) 
        ma <- ma[2]
    }
    if (is.matrix(x$var.summary[[i]])) 
      newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(x$var.summary[[i]]), 
                          byrow = TRUE)
    else newd[[i]] <- rep(ma, n.grid * n.grid)
  }
  names(newd) <- v.names
  newd[[view[1]]] <- v1
  newd[[view[2]]] <- v2
  if (type == "link") 
    zlab <- paste("linear predictor")
  else if (type == "response") 
    zlab <- type
  else stop("type must be \"link\" or \"response\"")
  fv <- predict.gam(x, newdata = newd, se.fit = TRUE, type = type)
  z <- fv$fit
  if (too.far > 0) {
    ex.tf <- exclude.too.far(v1, v2, x$model[, view[1]], 
                             x$model[, view[2]], dist = too.far)
    fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
  }
  if (is.factor(m1)) {
    m1 <- as.numeric(m1)
    m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
  }
  if (is.factor(m2)) {
    m2 <- as.numeric(m2)
    m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
  }
  if (se <= 0) {
    old.warn <- options(warn = -1)
    av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), n.grid, 
                 n.grid - 1)
    options(old.warn)
    max.z <- max(z, na.rm = TRUE)
    z[is.na(z)] <- max.z * 10000
    z <- matrix(z, n.grid, n.grid)
    surf.col <- t(av) %*% z %*% av
    surf.col[surf.col > max.z * 2] <- NA
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      min.z <- min(fv$fit, na.rm = TRUE)
      max.z <- max(fv$fit, na.rm = TRUE)
    }
    surf.col <- surf.col - min.z
    surf.col <- surf.col/(max.z - min.z)
    surf.col <- round(surf.col * nCol)
    con.col <- 1
    if (color == "heat") {
      pal <- heat.colors(nCol)
      con.col <- 3
    }
    else if (color == "topo") {
      pal <- topo.colors(nCol)
      con.col <- 2
    }
    else if (color == "cm") {
      pal <- cm.colors(nCol)
      con.col <- 1
    }
    else if (color == "terrain") {
      pal <- terrain.colors(nCol)
      con.col <- 2
    }
    else if (color == "gray" || color == "bw") {
      pal <- gray(seq(0.1, 0.9, length = nCol))
      con.col <- 1
    }
    else stop("color scheme not recognised")
    if (is.null(contour.col)) 
      contour.col <- con.col
    surf.col[surf.col < 1] <- 1
    surf.col[surf.col > nCol] <- nCol
    if (is.na(col)) 
      col <- pal[as.array(surf.col)]
    
    if(se.plot == TRUE)
      z.se <- matrix(fv$se.fit, n.grid, n.grid)
    else
      z.se <- matrix(fv$fit, n.grid, n.grid)
    z <- matrix(fv$fit, n.grid, n.grid)
    
    if (plot.type == "contour") {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",main=zlab"), ",...)", 
                    sep = "")
      if (color != "bw") {
        txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        
        eval(parse(text = txt))
        txt <- paste("contour(m1,m2,z.se,col=contour.col", 
                     ifelse("add" %in% dnm, "", ",add=TRUE"), ",...)", 
                     sep = "")
        eval(parse(text = txt))
      }
      else {
        txt <- paste("contour(m1,m2,z.se,col=1", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
    else {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",zlab=zlab"), ",...)", 
                    sep = "")
      if (color == "bw") {
        op <- par(bg = "white")
        txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", 
                     stub, sep = "")
        eval(parse(text = txt))
        par(op)
      }
      else {
        txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
  }
  else {
    if (color == "bw" || color == "gray") {
      subs <- paste("grey are +/-", se, "s.e.")
      lo.col <- "gray"
      hi.col <- "gray"
    }
    else {
      subs <- paste("red/green are +/-", se, "s.e.")
      lo.col <- "green"
      hi.col <- "red"
    }
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      z.max <- max(fv$fit + fv$se.fit * se, na.rm = TRUE)
      z.min <- min(fv$fit - fv$se.fit * se, na.rm = TRUE)
    }
    zlim <- c(z.min, z.max)
    z <- fv$fit - fv$se.fit * se
    z <- matrix(z, n.grid, n.grid)
    if (plot.type == "contour") 
      warning("sorry no option for contouring with errors: try plot.gam")
    stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                  ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
                                                                         dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, 
                                                                                                        "", ",sub=subs"), ",...)", sep = "")
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=lo.col"), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=\"black\""), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit + se * fv$se.fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=hi.col"), stub, sep = "")
    eval(parse(text = txt))
  }
}


#' @rdname diet-Internal
textG.dpart <- function (x, splits = TRUE, which = 4, label = "yval", FUN = text,
            all.leaves = FALSE, pretty = NULL, digits = getOption("digits") - 
              2, tadj = 0.65,  use.n = FALSE, bars = TRUE,
            xadj = 1, yadj = 1, bord = FALSE, node.cols = NULL, pos = NULL,
            ...) 
  {
    
    if (!inherits(x, "rpart")) 
      stop("Not legitimate rpart")
    #  if (!is.null(x$frame$splits)) 
    #    x <- rpconvert(x)  # commented this out as compiler didn't like this
    frame <- x$frame
    col <- names(frame)
    method <- x$method
    ylevels <- attr(x, "ylevels")
    if (!is.null(ylevels <- attr(x, "ylevels"))) 
      col <- c(col, ylevels)
    if (is.na(match(label, col))) 
      stop("Label must be a column label of the frame component of the tree")
    cxy <- par("cxy")
    if (!is.null(srt <- list(...)$srt) && srt == 90) 
      cxy <- rev(cxy)
    
    xy <- rpartco.dpart(x)
    node <- as.numeric(row.names(x$frame))
    is.left <- (node%%2 == 0)
    node.left <- node[is.left]
    parent <- match(node.left/2, node)
    
    if (splits) {
      left.child <- match(2 * node, node)
      right.child <- match(node * 2 + 1, node)
      rows <- labels(x, pretty = pretty)
      if (which == 1) 
        FUN(xy$x, xy$y + tadj * cxy[2], rows[left.child], 
            ...)
      else {
        if (which == 2 | which == 4) 
          FUN(xy$x, xy$y + tadj * cxy[2], rows[left.child], 
              pos = 2, ...)
        if (which == 3 | which == 4) 
          FUN(xy$x, xy$y + tadj * cxy[2], rows[right.child], 
              pos = 4, ...)
      }
    }
    leaves <- if (all.leaves) 
      rep(TRUE, nrow(frame))
    else frame$var == "<leaf>"
    
    if (is.null(frame$yval2))
      stat <- x$functions$text(yval = frame$yval[leaves],
                               dev = frame$dev[leaves], wt = frame$wt[leaves],
                               ylevel = ylevels, digits = digits, n = frame$n[leaves],
                               use.n = use.n)
    else stat <- x$functions$text(yval = frame$yval2[leaves,
                                                     ], dev = frame$dev[leaves], wt = frame$wt[leaves],
                                  ylevel = ylevels, digits = digits, n = frame$n[leaves],
                                  use.n = use.n)
    
    if(is.null(node.cols)){
      stat.cols <- as.factor(stat)
      levels(stat.cols) <- rainbow(length(levels(stat.cols)))
    }
    else{
      stat.cols <- as.factor(unlist(lapply(strsplit(stat, " "), function(x) x[1])))
      levels(stat.cols) <-  node.cols[match(levels(stat.cols), names(node.cols))]
    }
    stat <- as.factor(stat)
    
    # Plotting
    text.adj <- yadj * diff(range(xy$y))/50
    FUN(xy$x[leaves], xy$y[leaves] - tadj * cxy[2] - text.adj,
        stat, adj = 0.5, col = "black", ...)
    
    points(xy$x[leaves], xy$y[leaves], pch = 16, cex = 3 *
             par()$cex, col = as.vector(stat.cols))
    
    if(is.null(pos)){
      rx <- range(xy$x)
      ry <- range(xy$y)
      x.leg <- min(xy$x) - 0 * rx
      y.leg <- max(xy$y) + 0 * ry
    }
    
    if(is.null(node.cols)){
      if(is.null(pos))
        legend(x.leg, y.leg, levels(stat), col = levels(stat.cols),
               pch = 16, bty = "n", ...)
      else
        legend(pos, levels(stat), col = levels(stat.cols),
               pch = 16, bty = "n", ...)
    }
    else{
      if(is.null(pos))
        legend(x.leg, y.leg, names(node.cols), col = node.cols,
               pch = 16, bty = "n", ...)
      else
        legend(pos, names(node.cols), col = node.cols,
               pch = 16, bty = "n", ...)
    }
    
    invisible()
  }


#' @rdname diet-Internal
text.dpart <-
  function (x, splits = TRUE, which = 4, label = "yval", FUN = text,
            all.leaves = FALSE, pretty = NULL, digits = getOption("digits") - 
              2, tadj = 0.65,  use.n = FALSE, bars = TRUE,
            xadj = 1, yadj = 1, bord = FALSE, node.cols = NULL, pos = NULL,
            ...) 
  {
    
    if (!inherits(x, "rpart")) 
      stop("Not legitimate rpart")
    #   if (!is.null(x$frame$splits))  # commented out as compiler didn't like this
    #       x <- rpconvert(x)
    frame <- x$frame
    col <- names(frame)
    method <- x$method
    ylevels <- attr(x, "ylevels")
    if (!is.null(ylevels <- attr(x, "ylevels"))) 
      col <- c(col, ylevels)
    if (is.na(match(label, col))) 
      stop("Label must be a column label of the frame component of the tree")
    cxy <- par("cxy")
    if (!is.null(srt <- list(...)$srt) && srt == 90) 
      cxy <- rev(cxy)
    m <- match.call(expand.dots = TRUE)
    
    if(is.null(m$uniform))
      xy <- rpartco.dpart(x)
    else
      xy <- rpart:::rpartco(x)
    node <- as.numeric(row.names(x$frame))
    is.left <- (node%%2 == 0)
    node.left <- node[is.left]
    parent <- match(node.left/2, node)
    
    if (splits) {
      left.child <- match(2 * node, node)
      right.child <- match(node * 2 + 1, node)
      rows <- labels(x, pretty = pretty)
      if (which == 1) 
        FUN(xy$x, xy$y + tadj * cxy[2], rows[left.child], 
            ...)
      else {
        if (which == 2 | which == 4) 
          FUN(xy$x, xy$y + tadj * cxy[2], rows[left.child], 
              pos = 2, ...)
        if (which == 3 | which == 4) 
          FUN(xy$x, xy$y + tadj * cxy[2], rows[right.child], 
              pos = 4, ...)
      }
    }
    leaves <- if (all.leaves) 
      rep(TRUE, nrow(frame))
    else frame$var == "<leaf>"
    
    if (is.null(frame$yval2))
      stat <- x$functions$text(yval = frame$yval[leaves],
                               dev = frame$dev[leaves], wt = frame$wt[leaves],
                               ylevel = ylevels, digits = digits, n = frame$n[leaves],
                               use.n = use.n)
    else stat <- x$functions$text(yval = frame$yval2[leaves,
                                                     ], dev = frame$dev[leaves], wt = frame$wt[leaves],
                                  ylevel = ylevels, digits = digits, n = frame$n[leaves],
                                  use.n = use.n)
    
    if(is.null(node.cols)){
      stat.cols <- as.factor(stat)
      levels(stat.cols) <- rainbow(length(levels(stat.cols)))
    }
    else{
      stat.cols <- as.factor(unlist(lapply(strsplit(stat, " "), function(x) x[1])))
      levels(stat.cols) <-  node.cols[match(levels(stat.cols), names(node.cols))]
    }
    stat <- as.factor(stat)
    
    # Plotting
    text.adj <- yadj * diff(range(xy$y))/50
    FUN(xy$x[leaves], xy$y[leaves] - tadj * cxy[2] - text.adj,
        stat, adj = 0.5, col = "black", ...)
    
    points(xy$x[leaves], xy$y[leaves], pch = 16, cex = 3 *
             par()$cex, col = as.vector(stat.cols))
    
    if(is.null(pos)){
      rx <- range(xy$x)
      ry <- range(xy$y)
      x.leg <- min(xy$x) - 0 * rx
      y.leg <- max(xy$y) + 0 * ry
    }
    
    
    
    invisible()
  }


#' @rdname diet-Internal
tree.depth <-
  function (nodes)
  {
    depth <- floor(log(nodes, base = 2) + 1e-07)
    as.vector(depth - min(depth))
  }


#' @rdname diet-Internal
tree.depth.dpart <- function (nodes)
{
  depth <- floor(log(nodes, base = 2) + 1e-07)
  as.vector(depth - min(depth))
}


#' @rdname diet-Internal
snip.dpart <- function (x, toss)
{
  
  if (!inherits(x, "rpart")) 
    stop("Not an rpart object")
  if (missing(toss) || length(toss) == 0L) {
    toss <- snip.dpart.mouse(x)
    if (length(toss) == 0L) 
      return(x)
  }
  
  ff <- x$frame
  id <- as.numeric(row.names(ff))
  ff.n <- length(id)
  toss <- unique(toss)
  toss.idx <- match(toss, id, nomatch = 0)
  if (any(toss.idx == 0L)) {
    warning("Nodes ", toss[toss.idx == 0L], " are not in this tree")
    toss <- toss[toss.idx > 0L]
    toss.idx <- toss.idx[toss.idx > 0L]
  }
  id2 <- id
  while (any(id2 > 1)) {
    id2 <- floor(id2/2)
    xx <- (match(id2, toss, nomatch = 0) > 0)
    toss <- c(toss, id[xx])
    id2[xx] <- 0
  }
  temp <- match(floor(toss/2), toss, nomatch = 0)
  newleaf <- match(toss[temp == 0], id)
  keepit <- (1:ff.n)[is.na(match(id, toss))]
  n.split <- rep((1L:ff.n), ff$ncompete + ff$nsurrogate + 1 * 
                   (ff$var != "<leaf>"))
  split <- x$splits[match(n.split, keepit, nomatch = 0) > 0, 
                    , drop = FALSE]
  temp <- split[, 2L] > 1
  if (any(temp)) {
    x$csplit <- x$csplit[split[temp, 4L], , drop = FALSE]
    split[temp, 4] <- 1
    if (is.matrix(x$csplit)) 
      split[temp, 4L] <- 1L:nrow(x$csplit)
  }
  else x$csplit <- NULL
  x$splits <- split
  ff$ncompete[newleaf] <- ff$nsurrogate[newleaf] <- 0L
  ff$var[newleaf] <- "<leaf>"
  x$frame <- ff[sort(c(keepit, newleaf)), ]
  id2 <- id[x$where]
  id3 <- id[sort(c(keepit, newleaf))]
  temp <- match(id2, id3, nomatch = 0)
  while (any(temp == 0)) {
    id2[temp == 0] <- floor(id2[temp == 0]/2)
    temp <- match(id2, id3, nomatch = 0)
  }
  x$where <- match(id2, id3)
  
  x
}



#' @rdname diet-Internal
#' @import "grDevices"
#' @import "graphics"
snip.dpart.mouse <- function (tree, parms = paste(".rpart.parms", dev.cur(), sep = ".")) 
{
  
  xy <- rpartco.dpart(tree)
  toss <- NULL
  ff <- tree$frame
  if (exists(parms, envir = .GlobalEnv)) {
    parms <- get(parms, envir = .GlobalEnv)
    branch <- parms$branch
  }
  else branch <- 1
  node <- as.numeric(row.names(tree$frame))
  draw <- rpart.branch(xy$x, xy$y, node, branch)
  lastchoice <- 0
  while (length(choose <- identify(xy, n = 1, plot = FALSE)) > 
         0) {
    if (ff$var[choose] == "<leaf>") {
      cat("Terminal node -- try again\n")
      next
    }
    if (choose != lastchoice) {
      cat("node number:", node[choose], " n=", ff$n[choose], 
          "\n")
      cat("    response=", format(ff$yval[choose]))
      if (is.null(ff$yval2)) 
        cat("\n")
      else if (is.matrix(ff$yval2)) 
        cat(" (", format(ff$yval2[choose, ]), ")\n")
      else cat(" (", format(ff$yval2[choose]), ")\n")
      cat("    Error (dev) = ", format(ff$dev[choose]), 
          "\n")
      lastchoice <- choose
    }
    else {
      id <- node[choose]
      id2 <- node
      while (any(id2 > 1)) {
        id2 <- floor(id2/2)
        temp <- (match(id2, id, nomatch = 0) > 0)
        id <- c(id, node[temp])
        id2[temp] <- 0
      }
      temp <- match(id, node[ff$var != "<leaf>"], nomatch = 0)
      lines(c(draw$x[, temp]), c(draw$y[, temp]), col = 0L)
      toss <- c(toss, node[choose])
    }
  }
  toss
}


#' @rdname diet-Internal
#' @import "graphics"
select.tree.dpart <- function(object, se, nsplits){
    
    if(missing(se) && missing(nsplits))
      stop("Either standard error (se) or number of splits (nsplits) need to be defined for pruning.\n")
    else if(missing(se)){
      cptable <- data.frame(object$cptable)
      repeat{
        chk <- match(nsplits, cptable$nsplit)
        if(is.na(chk))
          nsplits <- nsplits + 1
        else
          break
      }
      cp <- with(cptable, CP[nsplit == round(nsplits)])
    }
    else{
      
      cptable <- data.frame(object$cptable)
      min.xerror <- min(cptable$xerror)
      min.xstd <- cptable$xstd[cptable$xerror == min.xerror][1]
      error <- min.xerror + min.xstd * se
      tree.se <- cptable$xerror <= error
      pick.tree <- (1:nrow(cptable))[tree.se][1]
      
      #     tt.cnt <- 0
      #      for(i in 1:length(tt)){
      #         if(tt[i] == TRUE)
      #            tt.cnt <- tt.cnt + 1
      #         else
      #            break
      #      }
      
      #      tt <- cptable$xerror < error
      #      tt.cnt <- (1:length(tt))[tt][1]
      
      if(pick.tree == 1){
        cp <- cptable$CP[1]
        cat(paste(se, "SE Tree has no splits.\n"))
      }
      else
        cp <- cptable$CP[pick.tree] 
      #    cp <- 0.5 * (cptable$CP[pick.tree] + cptable$CP[pick.tree - 1])
      
      
    }
    
    # produce a diagnostic plot showing the results from the cptable using the rsq.rpart function
    rn <- nrow(cptable[cptable$CP >= cp,])
    
    par(mfrow=c(2,1), mar = c(5, 4, 4, 1) + 0.1)
    rsq(object, rn = rn)
    par(mfrow=c(1,1), mar = c(5,4,4,2) + 0.1)
    
    cp
  }



#' @rdname diet-Internal
select.tree <- function(object, ...)
  UseMethod("select.tree")

#' @rdname diet-Internal
#' @import "graphics"
rsq.dpart <- function (x, rn) {
  if (!inherits(x, "dpart")) 
    stop("Not legitimate rpart")
  
  p.rpart <- printcp(x)
  xstd <- p.rpart[, 5L]
  xerror <- p.rpart[, 4L]
  rel.error <- p.rpart[, 3L]
  nsplit <- p.rpart[, 2L]
  method <- x$method
  plot(nsplit, 1 - rel.error, xlab = "Number of Splits", ylab = "R-square", 
       ylim = c(0, 1), type = "o")
  points(nsplit, 1 - xerror, type = "o", lty = 2)
  points(nsplit[rn], 1-xerror[rn], col = "blue", pch = 16)
  points(nsplit[rn], 1-rel.error[rn], col = "blue", pch = 16)
  legend("bottomright", c("Resubstitution", "X Relative", "Selected Tree"), lty = c(1:2, NA), pch = c(NA, NA, 16),
         col = c("black", "black", "blue"), bty = "n", cex = 0.8)
  ylim <- c(min(xerror - xstd) - 0.1, max(xerror + xstd) + 
              0.1)
  plot(nsplit, xerror, xlab = "Number of Splits", ylab = "X Relative Error", 
       ylim = ylim, type = "o")
  segments(nsplit, xerror - xstd, nsplit, xerror + xstd)
  points(nsplit[rn], xerror[rn], col = "blue", pch = 16)
  lines(c(-2, p.rpart[, 2L], max(p.rpart[, 2L]+2)), rep(min(p.rpart[, 4L]), 
                                                        nrow(p.rpart)+2), lty = 3, col = "grey")
  legend("topright", c("Minimum CV error", "Selected Tree"), lty = c(3, NA), pch = c(NA, 16), 
         col = c("grey", "blue"), bty = "n", cex = 0.8)
  
  
  invisible()
}

#' @rdname diet-Internal
rsq <- function(x, ...) UseMethod("rsq")

#' @rdname diet-Internal
#' @importFrom "grDevices" "dev.cur"
rpartco.dpart <- function (tree, parms = paste(".rpart.parms", dev.cur(), sep = "."))
  {
    
    frame <- tree$frame
    method <- tree$method
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    is.leaf <- (frame$var == "<leaf>")
    if (exists(parms, envir = .GlobalEnv)) {
      parms <- get(parms, envir = .GlobalEnv)
      uniform <- parms$uniform
      nspace <- parms$nspace
      minbranch <- parms$minbranch
    }
    else {
      uniform <- FALSE
      nspace <- -1
      minbranch <- 0.3
    }
    if (uniform) 
      y <- (1 + max(depth) - depth)/max(depth, 4)
    else {
      y <- dev <- frame$dev
      temp <- split(seq(node), depth)
      parent <- match(floor(node/2), node)
      sibling <- match(ifelse(node%%2, node - 1, node + 1), 
                       node)
      for (i in temp[-1]) {
        temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
        y[i] <- y[parent[i]] - temp2
      }
      fudge <- minbranch * diff(range(y))/max(depth)
      for (i in temp[-1]) {
        temp2 <- dev[parent[i]] - (dev[i] + dev[sibling[i]])
        haskids <- !(is.leaf[i] & is.leaf[sibling[i]])
        y[i] <- y[parent[i]] - ifelse(temp2 <= fudge & haskids, 
                                      fudge, temp2)
      }
      y <- y/(max(y))
    }
    x <- double(length(node))
    x[is.leaf] <- seq(sum(is.leaf))
    left.child <- match(node * 2, node)
    right.child <- match(node * 2 + 1, node)
    temp <- split(seq(node)[!is.leaf], depth[!is.leaf])
    for (i in rev(temp)) x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])
    if (nspace < 0) 
      return(list(x = x, y = y))
    compress <- function(me, depth) {
      lson <- me + 1
      x <- x
      if (is.leaf[lson]) 
        left <- list(left = x[lson], right = x[lson], depth = depth + 
                       1, sons = lson)
      else left <- compress(me + 1, depth + 1)
      rson <- me + 1 + length(left$sons)
      if (is.leaf[rson]) 
        right <- list(left = x[rson], right = x[rson], depth = depth + 
                        1, sons = rson)
      else right <- compress(rson, depth + 1)
      maxd <- max(left$depth, right$depth) - depth
      mind <- min(left$depth, right$depth) - depth
      slide <- min(right$left[1:mind] - left$right[1:mind]) - 
        1
      if (slide > 0) {
        x[right$sons] <- x[right$sons] - slide
        x[me] <- (x[right$sons[1]] + x[left$sons[1]])/2
        x <<- x
      }
      else slide <- 0
      if (left$depth > right$depth) {
        templ <- left$left
        tempr <- left$right
        tempr[1:mind] <- pmax(tempr[1:mind], right$right - 
                                slide)
      }
      else {
        templ <- right$left - slide
        tempr <- right$right - slide
        templ[1:mind] <- pmin(templ[1:mind], left$left)
      }
      list(left = c(x[me] - nspace * (x[me] - x[lson]), templ), 
           right = c(x[me] - nspace * (x[me] - x[rson]), tempr), 
           depth = maxd + depth, sons = c(me, left$sons, right$sons))
    }
    compress(1, 1)
    list(x = x, y = y)
  }


#' @rdname diet-Internal
plotG.dpart <- function(x, node.cols = NULL, pos = NULL, ...){
  
  
  #  plot(x, keep.margins = TRUE, ...)
  rpart:::plot.rpart(x, ...)
  #plot.rpart(x, ...)
  textG.dpart(x, xpd = NA, pretty = TRUE, splits = TRUE, node.cols = node.cols, 
              pos = pos,  ...)
  
  
  val <- rpartco.dpart(x)
  ff <- x$frame
  nodes <- as.numeric(row.names(ff))
  temp <- rpart.branch(val$x, val$y, nodes, branch=1)
  x1 <- c(temp$x)[seq(1, length(c(temp$x)), by = 5)]
  y1 <- c(temp$y)[seq(1, length(c(temp$y)), by = 5)]
  x2 <- c(temp$x)[seq(4, length(c(temp$x)), by = 5)]
  y2 <- c(temp$y)[seq(4, length(c(temp$y)), by = 5)]
  xvec <- c(x1, x2)
  yvec <- c(y1, y2)
  nodeLR <- c(temp$nodeL, temp$nodeR)
  
  id <- seq_along(xvec)
  nID <- nodeLR[id]
  rn <- row.names(ff)[-1]
  ff.var <- as.vector(ff$var[-1])[match(nID, rn)]
  
  text(xvec[id][ff.var == "<leaf>"], yvec[id][ff.var == "<leaf>"],
       nID[ff.var == "<leaf>"], cex = 0.7, col = "black")
  
  
  
}


#' @rdname diet-Internal
#' @importFrom "grDevices" "chull"
outside <- function(polyx, polyy, x, y)
  {
    # This function is used in Mask and locates any values inside the
    # convex hull.
    
    tt <- chull(c(x, polyx), c(y, polyy))
    return(any(tt == 1))
  }


#' @rdname diet-Internal
na.igam <- function (x) 
{
  Terms <- attr(x, "terms")
  if (!is.null(Terms)) 
    yvar <- attr(Terms, "response")
  else yvar <- 0L
  if (yvar == 0L) {
    xmiss <- is.na(x)
    keep <- (xmiss %*% rep(1, ncol(xmiss))) < ncol(xmiss)
  }
  else {
    xmiss <- is.na(x[-yvar])
    ymiss <- is.na(x[[yvar]])
    if (is.matrix(ymiss)) 
      keep <- ((xmiss %*% rep(1, ncol(xmiss))) < ncol(xmiss)) & 
      ((ymiss %*% rep(1, ncol(ymiss))) == 0)
    else keep <- ((xmiss %*% rep(1, ncol(xmiss))) < ncol(xmiss)) & 
      !ymiss
  }
  if (all(keep)) 
    x
  else {
    temp <- seq(keep)[!keep]
    names(temp) <- row.names(x)[!keep]
    class(temp) <- c("na.rpart", "omit")
    structure(x[keep, , drop = FALSE], na.action = temp)
  }
}


#' @rdname diet-Internal
na.dpart <- function (x) 
{
  Terms <- attr(x, "terms")
  if (!is.null(Terms)) 
    yvar <- attr(Terms, "response")
  else yvar <- 0L
  if (yvar == 0L) {
    xmiss <- is.na(x)
    keep <- (xmiss %*% rep(1, ncol(xmiss))) < ncol(xmiss)
  }
  else {
    xmiss <- is.na(x[-yvar])
    ymiss <- is.na(x[[yvar]])
    if (is.matrix(ymiss)) 
      keep <- ((xmiss %*% rep(1, ncol(xmiss))) < ncol(xmiss)) & 
      ((ymiss %*% rep(1, ncol(ymiss))) == 0)
    else keep <- ((xmiss %*% rep(1, ncol(xmiss))) < ncol(xmiss)) & 
      !ymiss
  }
  if (all(keep)) 
    x
  else {
    temp <- seq(keep)[!keep]
    names(temp) <- row.names(x)[!keep]
    class(temp) <- c("na.rpart", "omit")
    structure(x[keep, , drop = FALSE], na.action = temp)
  }
}

#' @rdname diet-Internal
#' @importFrom "grDevices" "chull"
mask <- function(gridx, gridy, x, y)
{
  
  # gridx = grid x-coordinate that provides the mask
  # gridy = grid y-coordinate that provides the mask
  # x = prediction of x-coordinate
  # y = prediction of y-coordinate
  
  
  hull <- chull(x, y)
  polyx <- x[hull]
  polyx <- mean(polyx) + 1.005 * (polyx - mean(polyx))
  polyy <- y[hull]
  polyy <- mean(polyy) + 1.005 * (polyy - mean(polyy))
  n <- length(gridx)
  tt <- numeric(n)
  for(i in 1:n)
    tt[i] <- outside(polyx, polyy, gridx[i], gridy[i])
  
  return(tt)
}



#' @rdname diet-Internal
#' @importFrom "stats" "reorder"
explore.bag <- function(object, node, cols = NULL, showtitle = FALSE, axis.side = 2, cex = 1.0, ylim){

if(!(axis.side == 2 | axis.side ==4))
  stop("Incorrect axis.side specified. Only 2 or 4 accepted.")
nc <- ncol(object$m)
if(missing(ylim))
  ylimit <- c(-0.05, 1.05)
else
  ylimit <- ylim
bID <- match(node, row.names(object$m))

bars <- barplot(object$m[bID,], plot = FALSE)
x <- as.vector(bars)


#nodevals <- data.frame(m = object$m[bID,], lci95 = object$lci95[bID,], 
#                       uci95 = object$uci95[bID,], prey = names(data.frame(object$m)))
nodevals <- data.frame(object$m[bID,], object$lci95[bID,], 
                       object$uci95[bID,], names(data.frame(object$m)))
names(nodevals) <- c("m", "lci95", "uci95", "prey")

tmp <- as.matrix(nodevals[,2:3])
id <- apply(tmp, 1, function(x){ all(x == 0)})
nodevals[id,1:3] <- NA

preyO <- 1:length(nodevals$prey)


p <- ggplot(nodevals, mapping = aes_string(x = reorder("prey", preyO), y = "lci95")) + 
  geom_segment(stat = "identity", aes_string(xend = "prey", yend = "uci95", 
                                      colour = reorder(cols, preyO)), 
               lineend = "butt", size = 1.5, 
               arrow = arrow(ends = "both", angle = 90, length = unit(0.1, "cm"),
                             type = "closed")) + 
  scale_color_manual(values = as.vector(cols), labels = names(cols), 
                     name = "Prey")  +
  ylim(ylim) + xlab("") + ylab("Bootstrapped Proportion") + 
  ggtitle(paste("Node", node)) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.margin = unit(c(1,1,1.5,1.2), "cm")) 

print(p)

data.frame(m = object$m[bID,], v = object$v[bID,], lci95 = object$lci95[bID,],
           uci95 = object$uci95[bID,])
}

#' @rdname diet-Internal
#' @importFrom "graphics" "par"
#' @importFrom "stats" "reorder"
explore <- function(object, pred, pred.where, loss = NULL, node, cols = NULL,
                    showtitle = FALSE, labels = TRUE, cex = 1.0, ylim){

  pred.node <- pred[pred.where == paste(node),]
  if(length(pred.node) == 0)
    print(cat(paste("No data at node", node, "\n")))
  else if(is.null(nrow(pred.node)))
    pred.node.m <- pred.node
  else
    pred.node.m <- apply(as.matrix(pred.node), 2, mean)
  
  if(showtitle)
    par(mar = c(6,4,4,2)+0.1)
  n <- length(pred.node.m)
  
  bars <- barplot(pred.node.m, plot = FALSE)
  x <- as.vector(bars)
  if(labels)
    nms <- names(pred.node.m)
  else
    nms <- NULL
  names(pred.node.m) <- NULL
  if(missing(cols))
    cols <- "grey"
  
  if(missing(ylim))
    ylim <- c(-0.05,1.05)
  

  preyO <- 1:length(pred.node.m)
  p <- ggplot(mapping = aes(x = reorder(names(data.frame(pred.node)),preyO), 
                            y = pred.node.m)) + 
    geom_bar(stat = "identity", aes(fill = reorder(cols, preyO))) + 
    scale_fill_manual(values = as.vector(cols), labels = names(cols), name = "Prey")  +
    ylim(ylim) + xlab("") + ylab("Proportion") + 
    ggtitle(paste("Node ", row.names(object$frame)[as.integer(node)], "\n",
                  "Diet Composition (D=", round(loss, 3), ")", sep = "")) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          plot.margin = unit(c(1,1,1.5,1.2), "cm")) 
  
  print(p)
  
  
}



#' @rdname diet-Internal
#' @importFrom "utils" "write.csv"
writepn.csv <- function(x){
  
  if (!inherits(x, "diet"))
    stop("Not an object of class diet")
  
  lvls <- length(levels(x$Group))
  dat <- data.frame(Group = sort(levels(x$Group)), PreyTaxonSort = rep("", lvls), 
                    PreyTaxBroad = rep("", lvls), Sort = rep("", lvls))
  write.csv(dat, "PreyTaxonSort_Template.csv", row.names = FALSE)
  cat("Prey groupings written to file: PreyTaxonSort.csv, Please edit file to assign prey
      groupings. \n")
  
  
}



#' @rdname diet-Internal
#' @importFrom "mgcv" "gam"
#' @import "sp"
#' @import "spaMM" 
#' @importFrom "raster" "raster"
#' @import "lattice"
#' @importFrom "rasterVis" "levelplot"
#' @importFrom "utils" "data"
SmPlots <- function(x, i, SmXvar, SmXdat, LonID, LatID, projection, palette, too.far = 0.05){
  
  tmpdat <- data.frame(SmXdat[,i], x[,LonID], x[,LatID])
  names(tmpdat) <- c("y", "Longitude", "Latitude")
  
  #tmpdat <- data.frame(y = SmXdat[,i], Longitude = x[,LonID], Latitude = x[,LatID])
  res <- gam(y ~ s(tmpdat$Longitude, tmpdat$Latitude), data = tmpdat)

 # data(worldcountries)
  plotfit <- Vis.Gam(res, too.far = too.far, plot.type = "contour", add = FALSE)
  mat <- plotfit$mat
  names(mat) <- c("Longitude", "Latitude", "z")
  
  coordinates(mat) <- ~Longitude+Latitude  # coordinates are being set for the raster
  proj4string(mat) <- CRS(projection)  # projection is being set for the raster
  gridded(mat) <- TRUE  # a gridded structure is being set for the raster
  mat.raster <- raster(mat)  # the raster is being created
  
  zext <- range(mat$z, na.rm = TRUE)
  # Set up of plot

  p <- levelplot(mat.raster, maxpixels=4e6, margin=FALSE, cuts=length(palette)-1, col.regions=palette,
                 xlab = "", ylab = "", at = seq(zext[1], zext[2], length = 50),
                 main = SmXvar[i])
  
  # country layer
  country.layer <- layer(
    sp.polygons(worldcountries, fill=fill, col = col),
    data=list(sp.polygons=sp.polygons, worldcountries=worldcountries, 
              fill="darkgray", col = "lightgray") 
  )
  
  # points layer
  points.layer <- layer(
    panel.points(Longitude, Latitude, pch = 16, col = "black", cex = 0.6),
    data = tmpdat
  )
  
  
  smplot <- p + country.layer + points.layer
  
  list(res = res, smplot = smplot)
}


#' @rdname diet-Internal
#' @import "geoR"
Distance <- function(O, P, type = "Hellinger"){
  
  # O = observed matrix
  # P = predicted matrix
  
  W <- rowSums(O)
  rW <- W/sum(W)
  
  
  if(type == "KL"){
    KLDist <- function(O,E) rowSums(O*log((O + (O==0))/E))
    d <- sum(rW * (KLr <- KLDist(O, P)))
    val <- list(Dist = KLr, d = d)
    
  }
  else if(type == "Hellinger"){
    HDist <- function(O, E) sqrt(0.5*rowSums((sqrt(O) - sqrt(E))^2))
    d <- sum(rW * (Hr <- HDist(O, P)))
    val <- list(Dist = Hr, d = d)
  }
  
  val
}

