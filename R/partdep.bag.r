partdep.bag <- function (object, Xvar, Yvar = NULL, fact = FALSE, var.cond = NULL, 
    plotmap = FALSE, mfrow = NULL, xlim = NULL, ylim = NULL, 
    leg.pos = "topleft", plot2file = FALSE, se.fit = FALSE, too.far = 0.1,
    sp.id = NULL) 

{
    if (plot2file) 
        pdf("partdep.pdf")
    def.par <- par(no.readonly = TRUE)
    
    par(new = FALSE, mar = c(5, 4, 4, 2) + 0.1, mfrow = c(1, 
        1))
    if (is.null(mfrow)) 
        mfrow <- c(3, 3)
    if (!inherits(object, "bag")) 
        stop("Not a bagged object")
    if (length(var.cond) > 1) 
        stop("Only one conditioning variable can be specified.\n")
    dat <- object$data

    if (!is.null(var.cond)) {
        m <- match(var.cond[[1]], na.omit(unique(dat[, names(var.cond)])))
        if (length(m) == length(na.omit(unique(dat[, names(var.cond)])))) {
            var.cond <- NULL
            warning("Defaulting to var.cond=NULL as all categories in the conditioning variable\n         have been specified\n")
        }
    }
    if (length(Xvar) > 1) {
        if (missing(var.cond) | is.null(names(var.cond))) {
            warning("Defaulting to selecting the mean or median for each conditioning term as\n            var.cond is not specified. See ?partdep for more information.")
            Xcond <- NULL
            dat.cond <- dat
            x.labels <- attr(object$baggs[[1]]$terms, "term.labels")
            Xmat <- dat[, x.labels]
            Xmat <- Xmat[, -match(Xvar, names(Xmat))]
            Xcond <- data.frame(matrix(NA, ncol = ncol(Xmat), 
                nrow = 20 * 20))
            names(Xcond) <- names(Xmat)
            
            for (i in 1:ncol(Xmat)) {
                colX <- Xmat[, i]
                if (is.factor(colX)) {
                  Xcond[, i] <- mean(as.numeric(colX), na.rm = TRUE)
                  Xcond[, i] <- as.factor(Xcond[, i])
                  levels(Xcond[, i]) <- levels(colX)
                }
                else Xcond[, i] <- mean(colX, na.rm = TRUE)
            }
            X.id <- match(names(Xcond), x.labels)
            list.val <- lapply(Xcond, length)
            list.len <- max(unlist(list.val))
            nseq <- 20
            seqval <- NULL
            for (i in 1:nseq) seqval <- c(seqval, rep(i, nseq))
            if (is.factor(object$data[, match(Xvar[1], names(object$data))])) {
                Xseq <- as.vector(unique(object$data[, match(Xvar[1], 
                  names(object$data))]))
                Yseq <- as.vector(unique(object$data[, match(Xvar[2], 
                  names(object$data))]))
            }
            else {
                Xseq <- seq(min(object$data[, match(Xvar[1], 
                  names(object$data))], na.rm = TRUE), max(object$data[, 
                  match(Xvar[1], names(object$data))], na.rm = TRUE), 
                  length = nseq)
                Yseq <- seq(min(object$data[, match(Xvar[2], 
                  names(object$data))], na.rm = TRUE), max(object$data[, 
                  match(Xvar[2], names(object$data))], na.rm = TRUE), 
                  length = nseq)
                newdatX <- rep(Xseq, nseq)
                newdatY <- Yseq[seqval]
                newdat <- data.frame(cbind(newdatX, newdatY, 
                  Xcond))
                names(newdat) <- c(paste(Xvar), names(Xcond))
                pred <- predict(object$baggs[[1]], newdat, type = "prob", 
                  plot = FALSE)
                for (i in 2:length(object$baggs)) pred <- pred + 
                  predict(object$baggs[[i]], newdat, type = "prob", 
                    plot = FALSE)
                predc <- pred/length(object$baggs)
                if (!is.null(Yvar)) {
                  if (all(is.integer(Yvar)) | length(na.omit(match(Yvar, 
                    names(predc)))) > 0) 
                    predc <- predc[, Yvar]
                  else warning("Yvar not matched. Plots will be produced for all Yvars.\n")
                }
            }
        }
        else {

            x.labels <- attr(object$baggs[[1]]$terms, "term.labels")
            Xmat <- dat[, x.labels]
            Xmatn <- Xmat[, -c(match(Xvar[1], names(Xmat)), match(Xvar[2], 
                names(Xmat)))]
            tmp <- Xmatn[, match(names(var.cond), names(Xmatn))]
            id <- match(tmp, var.cond[[names(var.cond)]])
            Xmatn <- Xmatn[id & !is.na(id), ]
            dat.cond <- dat[id & !is.na(id), ]
            id <- match(names(var.cond), names(Xmatn))
            Xcondv <- data.frame(matrix(NA, ncol = ncol(Xmatn) - 
                length(id), nrow = length(var.cond[[1]])))
            names(Xcondv) <- names(Xmatn)[-id]
            row.names(Xcondv) <- var.cond[[1]]
            for (i in 1:ncol(Xmatn[, -id])) {
                colX <- Xmatn[, -id][, i]
                if (is.factor(colX)) {
                  val <- tapply(as.vector(colX), Xmatn[, id], 
                                function(y) mean(as.numeric(y), 
                                  na.rm = TRUE))
                  val <- val[!is.na(val)]
                  m <- match(names(val), row.names(Xcondv))
                  Xcondv[m, i] <- val
                  Xcondv[, i] <- as.factor(Xcondv[, i])
                  levels(Xcondv[, i]) <- levels(colX)
                }
                else {
                  val <- tapply(colX, as.vector(Xmatn[, id]), 
                    function(y) mean(y, na.rm = TRUE))
                  m <- match(names(val), row.names(Xcondv))
                  Xcondv[m, i] <- val
                }
            }
            if (length(var.cond[[names(var.cond)]]) == 1) {
                tmp <- data.frame(matrix(c(Xcondv, var.cond[[names(var.cond)]]), 
                  nrow = 1, byrow = TRUE))
                names(tmp) <- c(names(Xcondv), names(var.cond))
                Xcondv <- tmp
            }
            else {
                Xcondv <- data.frame(cbind(Xcondv, var.cond[[names(var.cond)]]))
                names(Xcondv)[ncol(Xcondv)] <- names(var.cond)
            }
            nseq <- 20
            seqval <- NULL
            for (i in 1:nseq) seqval <- c(seqval, rep(i, nseq))
            Xseq <- seq(min(object$data[, match(Xvar[1], names(object$data))], 
                na.rm = TRUE), max(object$data[, match(Xvar[1], 
                names(object$data))], na.rm = TRUE), length = nseq)
            Yseq <- seq(min(object$data[, match(Xvar[2], names(object$data))], 
                na.rm = TRUE), max(object$data[, match(Xvar[2], 
                names(object$data))], na.rm = TRUE), length = nseq)
            newdatX <- rep(Xseq, nseq)
            newdatY <- Yseq[seqval]
            if (nrow(Xcondv) > 1) {
                predc <- list()
                for (j in 1:nrow(Xcondv)) {
                  if (!any(is.na(Xcondv[j, ]))) {
                    newdat <- data.frame(cbind(newdatX, newdatY, 
                      matrix(rep(Xcondv[j, ], length(newdatX)), 
                        byrow = TRUE, nrow = length(newdatX))))
                    names(newdat) <- c(paste(Xvar), names(Xcondv))
                    for (i in 1:ncol(newdat)) {
                      m <- match(names(newdat)[i], names(object$data))
                      if (is.logical(object$data[, m])) 
                        newdat[, i] <- as.logical(newdat[, i])
                      else if (is.numeric(object$data[, m])) 
                        newdat[, i] <- as.numeric(newdat[, i])
                      else if (is.ordered(object$data[, m])) 
                        newdat[, i] <- as.ordered(newdat[, i])
                      else if (is.factor(object$data[, m])) 
                        newdat[, i] <- as.factor(unlist(newdat[, 
                          i]))
                    }
                    pred <- predict(object$baggs[[1]], newdat, 
                      type = "prob", plot = FALSE)
                    for (i in 2:length(object$baggs)) pred <- pred + 
                      predict(object$baggs[[i]], newdat, type = "prob", 
                        plot = FALSE)
                    predc[[j]] <- pred/length(object$baggs)
                  }
                  else predc[[j]] <- NULL
                }
                sumpredc <- matrix(0, nrow = nrow(newdat), ncol = ncol(object$pred[[1]]))
                for (i in 1:length(predc)) if (!is.null(predc[[i]])) 
                  sumpredc <- predc[[i]] + sumpredc
                predc <- sumpredc/length(predc[!sapply(predc, 
                  is.null)])
            }
            else {
                tmp <- unlist(Xcondv)
                newdat <- data.frame(cbind(newdatX, newdatY, 
                  matrix(rep(tmp, length(newdatX)), byrow = TRUE, 
                    nrow = length(newdatX))), stringsAsFactors = FALSE)
                names(newdat) <- c(paste(Xvar), names(Xcondv))
                cl <- attr(delete.response(object$baggs[[1]]$terms), 
                  "dataClasses")
                id <- match(names(newdat), names(cl))
                for (i in 1:length(cl[id])) {
                  if (cl[id][i] == "numeric") 
                    newdat[, names(cl[id])[i]] <- as.numeric(newdat[, 
                      names(cl[id])[i]])
                  else if (cl[id][i] == "factor") 
                    newdat[, names(cl[id])[i]] <- as.factor(newdat[, 
                      names(cl[id])[i]])
                  else if (cl[id][i] == "logical") 
                    newdat[, names(cl[id])[i]] <- as.logical(newdat[, 
                      names(cl[id][i])])
                }
                pred <- predict(object$baggs[[1]], newdat, type = "prob", 
                  plot = FALSE)
                for (i in 2:length(object$baggs)) pred <- pred + 
                  predict(object$baggs[[i]], newdat, type = "prob", 
                    plot = FALSE)
                predc <- pred/length(object$baggs)

            }
            if (!is.null(Yvar)) {
                if (all(is.integer(Yvar)) | length(na.omit(match(Yvar, 
                  names(predc)))) > 0) 
                  predc <- predc[, Yvar]
                else warning("Yvar not matched. Plots will be produced for all Yvars.\n")
            }
        }

       x <- as.numeric(as.vector(newdat[, Xvar[1]]))
        y <- as.numeric(as.vector(newdat[, Xvar[2]]))
        varx <- dat[, Xvar[1]]
        vary <- dat[, Xvar[2]]
        if(!is.null(var.cond)){
         varx.pts <- varx[dat[,names(var.cond)] %in% var.cond[[1]]]
          vary.pts <- vary[dat[,names(var.cond)] %in% var.cond[[1]]]
        }
        else{
          varx.pts <- varx
          vary.pts <- vary
        }
      if (too.far > 0) 
          tt <- exclude.too.far(x, y, varx, vary, dist = too.far)

        tt <- ifelse(tt == 1, TRUE, FALSE)
        nprey <- ncol(predc)
        value <- matrix(predc[, sp.id], nrow = nseq, 
                  ncol = nseq)
                value[tt] <- NA
                if(plotmap){
                  
                  projection <- "+proj=longlat +datum=WGS84"
                  
                  mat <- expand.grid(Longitude = Xseq, Latitude = Yseq)
                  mat$p <- as.vector(value)
                  
                  
                  
                  coordinates(mat) <- ~Longitude+Latitude  # coordinates are being set for the raster
                  proj4string(mat) <- CRS(projection)  # projection is being set for the raster
                  gridded(mat) <- TRUE  # a gridded structure is being set for the raster
                  mat.raster <- raster(mat)  # the raster is being created
                  
                  # Set up of plot
                  p <- levelplot(mat.raster, maxpixels=4e6, margin=FALSE, cuts=length(palette)-1, col.regions=pal,
                                 xlab = "", ylab = "", at = seq(0, 1, length = 50))
                  
                  # country layer
                  country.layer <- layer(
                    sp.polygons(worldcountries, fill=fill, col = col),
                    data=list(sp.polygons=sp.polygons, worldcountries=worldcountries, 
                              fill="lightgrey", col = "darkgrey") 
                  )
                  
                  # points layer
                  dat.pts <- list(x = varx.pts, y = vary.pts)
                  points.layer <- layer(
                    panel.points(x, y, ch = 16, col = "black", cex = 0.6),
                    data = dat.pts
                  )
                  
                  
                  print(p + country.layer + points.layer)
               }
    }
    else {
        pred <- predict(object$baggs[[1]], dat, type = "prob", 
            plot = FALSE)
        for (i in 2:length(object$baggs)) pred <- pred + predict(object$baggs[[i]], 
            dat, type = "prob", plot = FALSE)
        predc <- pred/length(object$baggs)
        X <- dat[, Xvar]
        if (fact) 
            X <- as.factor(X)
        if (length(Xvar) > 1) 
            predcN <- predc[!apply(X, 1, function(x) any(is.na(x))), 
                ]
        else predcN <- predc[!is.na(X), ]
        X <- na.omit(X)
        options(warn = 0)
        if (!is.null(Yvar)) {
            if (all(is.integer(Yvar)) | length(na.omit(match(Yvar, 
                names(predcN)))) > 0) 
                predcN <- predcN[, Yvar]
            else warning("Yvar not matched. Plots will be produced for all Yvars.\n")
        }
        if (is.null(ncol(predcN))) 
            predcN <- matrix(predcN, byrow = TRUE, ncol = 1)
        mat <- mat.L <- mat.U <- data.frame(matrix(NA, ncol = ncol(predcN), 
            nrow = length(unique(X))))
        for (i in 1:ncol(predcN)) {
            mat[, i] <- tapply(predcN[, i], X, mean)
            mat.L[, i] <- tapply(predcN[, i], X, function(x) quantile(x, 
                0.1))
            mat.U[, i] <- tapply(predcN[, i], X, function(x) quantile(x, 
                0.9))
        }
        names(mat) <- names(predcN)
        names(mat.L) <- names(predcN)
        names(mat.U) <- names(predcN)
        if (is.null(ylim)) 
            ylim <- c(0, 1)
        if (is.factor(X)) {
            nms <- names(mat)
            par(mfrow = mfrow, mar = c(5, 4, 4, 2) + 0.1, ask = TRUE)
            if (se.fit) {
                options(warn = -1)
                for (i in 1:ncol(mat)) {
                  plot(1:nrow(mat), mat[, i], pch = 15, ylim = c(0, 
                    1), main = nms[i], ylab = "Proportion", axes = FALSE, 
                    xlab = "")
                  axis(side = 2)
                  box()
                  arrows(1:nrow(mat), mat.L[, i], 1:nrow(mat), 
                    mat.U[, i], length = 0.1, angle = 90, code = 3)
                  nmsX <- levels(X)
                  if (any(sapply(nmsX, nchar) > 3)) 
                    mtext(side = 1, nmsX, at = 1:nrow(mat), adj = 1, 
                      cex = 0.6, las = 2, line = 1)
                  else mtext(side = 1, nmsX, at = 1:nrow(mat), 
                    adj = 0.5, cex = 0.6, , line = 1)
                }
                options(warn = 0)
            }
            else {
                for (i in 1:ncol(mat)) {
                  val <- barplot(mat[, i], names.arg = "", main = nms[i], 
                    ylim = ylim, ylab = "Proportion", axes = FALSE)
                  axis(side = 2)
                  x.nms <- levels(X)
                  if (any(sapply(x.nms, nchar) > 3)) 
                    mtext(side = 1, x.nms, at = val, adj = 1, 
                      cex = 0.6, las = 2, line = 1)
                  else mtext(side = 1, x.nms, at = val, adj = 0.5, 
                    cex = 0.6, line = 1)
                }
            }
        }
        else {
            if (is.integer(Xvar)) 
                xlabel <- names(newdat)[Xvar]
            else xlabel <- Xvar
            lun <- length(unique(X))
            xseq <- seq(min(X), max(X), length = ifelse(lun > 
                10, 10, 5))
            Xnew <- tapply(X, cut(X, xseq, include.lowest = TRUE), 
                mean)
            matN <- matL <- matU <- data.frame(matrix(NA, ncol = ncol(predcN), 
                nrow = length(Xnew)))
            names(matN) <- names(predcN)
            names(matL) <- names(predcN)
            names(matU) <- names(predcN)
            for (i in 1:ncol(predcN)) {
                matN[, i] <- tapply(predcN[, i], cut(X, xseq, 
                  include.lowest = TRUE), mean)
                matL[, i] <- tapply(predcN[, i], cut(X, xseq, 
                  include.lowest = TRUE), function(x) quantile(x, 
                  0.1))
                matU[, i] <- tapply(predcN[, i], cut(X, xseq, 
                  include.lowest = TRUE), function(x) quantile(x, 
                  0.9))
            }
            nms <- names(mat)
            matN <- na.omit(matN)
            matL <- na.omit(matL)
            matU <- na.omit(matU)
            Xnew <- Xnew[!is.na(Xnew)]
            par(mfrow = mfrow, mar = c(5, 4, 4, 2) + 0.1, ask = TRUE)
            options(warn = -1)
            for (i in 1:ncol(mat)) {
                plot(Xnew, matN[, i], type = "l", main = nms[i], 
                  ylim = ylim, xlab = xlabel, ylab = "Proportion")
                if (se.fit) {
                  xx <- c(Xnew, rev(Xnew))
                  yy <- c(matL[, i], rev(matU[, i]))
                  polygon(xx, yy, col = "light grey", border = NA)
                  lines(Xnew, matN[, i])
                }
                rug(X[dat$Group == levels(dat$Group)[i]])
            }
        }
        options(warn = 0)
        par(mfrow = c(1, 1), ask = FALSE, mar = c(5, 4, 4, 2) + 
            0.1)
    }
    if (plot2file) 
        dev.off()
    
    par(def.par)
}
