predict.dpart <-
function(object, newdata = list(), type = c("prob", "class"), na.action = na.pass,
   plot = TRUE, pred.type = "obs", predatorID = NULL, cex = 1.0, ...){

   if(plot)
      def.par <- par(no.readonly = TRUE)
  
   if (!inherits(object, "dpart"))
        stop("Not legitimate tree")
   if(!(type == "prob" | type == "class"))
     stop("Incorrect type specified. Either 'prob' or 'class' required.")
        
   if(pred.type == "predator" & is.null(predatorID))
      stop("Predator ID not specified.")

   if (missing(newdata)){
        where <- object$where
        ID <- object$data[,paste(predatorID)]
   }
   else{
       if (is.null(attr(newdata, "terms"))) {
          Terms <- delete.response(object$terms)
          ID <- newdata[,predatorID]
          newdata <- model.frame(Terms, newdata, na.action = na.action, 
                              xlev = attr(object, "xlevels"))
          if (!is.null(cl <- attr(Terms, "dataClasses"))) 
            .checkMFClasses(cl, newdata, TRUE)
      }
     else
       ID <- newdata[,paste(predatorID)]
    where <- rpart:::pred.rpart(object, rpart.matrix(newdata))
        
   }

   frame <- object$frame
   method <- object$method
   ylevels <- attr(object, "ylevels")
   nclass <- length(ylevels)
  
   if (type == "class" && nclass > 0) {
       pred <- factor(ylevels[frame$yval[where]], levels = ylevels)
       names(pred) <- names(where)
       tpred <- table(pred)/sum(table(pred))
       if(plot){
       
          if(dev.cur() == 1) windows(8,8, record = TRUE)
             par.old <- par()$mar
          par(mfrow=c(1,1), mar = c(8, 4, 4, 2) + 0.1)
          b <- barplot(tpred, names = "",
                main = "Distribution of predicted classifications at terminal nodes")
                mtext(side = 1, text = levels(pred), at = as.vector(b), las = 2, adj = 1, line = 1)
          par(mar = par.old)
       }
   }
   else if (type == "prob" && nclass > 0L) {
             pred <- frame$yval2[where, 1L + nclass + 1L:nclass, drop = FALSE]
	           dimnames(pred) <- list(names(where), ylevels)
	           if(plot){
                if(dev.cur() == 1) windows(8,8, record = TRUE)
                par(mfrow=c(2,2))
                for(i in names(table(where)))
                    explore(object = object, pred=pred, pred.where=where,
                       node=i, showtitle = TRUE, cex = cex)
                par(mfrow=c(1,1))
             }
          }

   pred <- data.frame(pred)
   
   if(pred.type == "predator"){
      pred$ID <- ID
      pred$where <- as.integer(where)
      if(type == "class")
         pred <- apply(pred, 2, function(x, ID) tapply(x, ID, unique),  pred$ID)
      else
         pred <- apply(pred, 2, function(x, ID) tapply(x, ID, mean),  pred$ID)
      pred <- data.frame(pred)
      pred <- pred[,-match(c("ID", "where"), names(pred))]
   }
   if(plot)
      par(def.par)
  
   pred

}




