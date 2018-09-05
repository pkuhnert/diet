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