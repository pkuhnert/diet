print.dpart <- function(x, setID = NULL, digits=getOption("digits"), file = "diet_tree_summary.csv", ...){

  object <- x
   ff <- object$frame
   node <- as.integer(row.names(ff))
   nobs <- ff$n
   npredators <- ff$wt
   dev <- ff$dev
   expLoss <- ff$dev/npredators

   val <- rpartco.dpart(object)
   temp <- rpart.branch(val$x, val$y, node, branch=1)
   nodeLR <- sort(c(1, temp$nodeL, temp$nodeR))
   nprey <- nsets <- pclass <- splits <- NULL
   for(i in 1:length(nodeLR)){
     subtree <- snip.dpart(object, nodeLR[i])
     nodesub <- as.numeric(row.names(subtree$frame))
     val <- rpartco.dpart(subtree)
     rn <- row.names(subtree$frame)[subtree$frame$var == "<leaf>"]
     pred.where <- names(table(subtree$where))[(1:length(rn))[rn == nodeLR[i]]]
     subpred <- predict(subtree, type = "prob", plot = FALSE)
     subpred.where <- subtree$where
     pred.node <- subpred[subpred.where == paste(pred.where),]
     pred.node.m <- apply(pred.node, 2, mean)
     nprey[i] <- length(pred.node.m[pred.node.m > 0])
     if(is.null(setID))
        nsets[i] <- NA
     else
        nsets[i]  <- length(unique(subtree$data[subpred.where ==
           as.integer(pred.where),][,setID]))
    }
    # Extracting the predicted class
    pclass <- levels(object$data$Group)[ff$yval]
    
    # Extracting the split information
    splits <- ff$var[order(node)]
    options(warn = -1) # suppresses an unnecessary warning message

    xc <- rpconvert(object)
   # xc <- object
    options(warn = 0)
    xc.sp <- xc$splits

   cuts <- vector(mode='character', length=nrow(xc$splits))
   temp <- xc$splits[ ,2L]
   for (i in 1L:length(cuts)) {
     if(temp[i] !=0){ 
       if (temp[i] == -1L)
         cuts[i] <-paste("<", format(signif(xc$splits[i,4], digits=digits)))
       else if (temp[i] ==1L)
         cuts[i] <-paste(">", format(signif(xc$splits[i,4], digits=digits)))
       else cuts[i]<- paste("splits as ",
                            paste(c("L", "-", "R")[xc$csplit[xc$splits[i,4], 1L:temp[i]]],
                                  collapse='', sep=''), collapse='')
     }
   }
    
   
    splitval <- vector("numeric", length = nrow(ff))
    for(i in 1:nrow(ff)){    
        m <- (row.names(xc$splits) == ff$var[i]) & 
          (xc$splits[,1] == ff$n[i]) 
        if(length(cuts[m]) !=0){
          splitval[i] <- cuts[m][1]
        }
        else
          splitval[i] <- NA
    }
    splitval <- splitval[order(node)]
    splitcomb <- paste(splits, splitval)
    splitcomb[is.na(splitval)] <- "<leaf>"

   sumdata <- data.frame(node = sort(node), nobs = nobs[order(node)], nsets = nsets,
         npredators = npredators[order(node)], nprey = nprey, dev = dev[order(node)], 
         loss = expLoss[order(node)], split = splitcomb, pclass = pclass[order(node)])

   write.csv(sumdata, file, row.names = FALSE)
   cat(paste("Written summary table out to file: ", file, "\n", sep = ""))

   invisible(sumdata)

}



