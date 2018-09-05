partdepS.bag <- function (object, Xvar, Yvar = NULL, fact = FALSE, var.cond = NULL, 
    plotmap = FALSE, mfrow = NULL, xlim = NULL, ylim = NULL, 
    leg.pos = "topleft", plot2file = FALSE, se.fit = FALSE, too.far = 0.1,
    sp.id = NULL) 
{
  
  data(worldcountries)  # loading from package maps
  pal <- spaMM.colors()
  
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
    if(is.null(sp.id))
      stop("Prey species not specified.")
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
                    function(y) mean(as.numeric(y), na.rm = TRUE))
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
        if (!is.null(var.cond)) {
            varx.pts <- varx[dat[, names(var.cond)] %in% var.cond[[1]]]
            vary.pts <- vary[dat[, names(var.cond)] %in% var.cond[[1]]]
        }
        else {
            varx.pts <- varx
            vary.pts <- vary
        }
      
       # tt <- mask(x, y, varx, vary)
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

    if (plot2file) 
        dev.off()
    par(def.par)
}

