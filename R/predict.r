#' Predict from a dpart Object
#' 
#' @description Forms predictions from an object of class \code{dpart}. 
#' Options are to produce predicted
#' probabilities or predicted classificaitons.
#' 
#' 
#' @param object object of class \code{dpart}
#' @param newdata optional data set to predict on.
#' @param type either "prob" or "class" to specify probabilities or 
#' classifications.
#' @param node.cols node.cols
#' @param plot logical. Produce a plot of the predictions if TRUE 
#' otherwise just write the predictions
#' out to file.
#' @param pred.type either "obs" or "predator" to indicate whether 
#' predictions are to be generated for
#' each observation or each predator, respectively.
#' @param predatorID predator identification label. Required when 
#' \code{pred.type = "predator"}
#' @param cex numeric. Size of plotting characters and labels.
#' @param \dots other arguments to pass into the function
#' 
#'                             
#'    
#' @return Predicted probabilities for each observation or 
#' each predator, if \code{pred.type = "predator"} was specified. 
#' Predicted classifications if \code{type = "class"} was specified.
#' 
#' @references Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The 
#' Analysis of Predator Diet and Stable Isotope Data, Journal of Statistical 
#' Software, In Prep. 
#' 
#' @examples 
#' # Load data
#' #data(yftdiet)  
#' # Load the prey taxa data
#' #data(PreyTaxonSort)
#' 
#' # Assigning prey colours for default palette
#' #val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = TRUE)
#' #node.colsY <- val$cols
#' #dietPP <- val$x   # updated diet matrix with Group assigned prey 
#' #taxa codes
#' 
#' # Fitting the classification tree
#' #yft.dp <- dpart(Group ~ Lat + Lon + Year + Quarter + SST  + Length,
#' #                   data = dietPP, weights = W, minsplit = 10,
#' #                                     cp = 0.001)
#'                                     
#' # Pruning the tree
#' #yft.pr <- prune(yft.dp, se = 1)                   
#' 
#' # Predictions
#' # predict distribution of prey composition for each predator
#' #yft.predator <- predict(yft.pr, type = "prob", pred.type = "predator",
#' #                          predatorID = "TripSetNo")
#'                           
#' # predict distribution of prey composition for each observation
#' #yft.pred.obs <- predict(yft.pr, type = "prob")
#'                           
#' # predict classification  for each observation in the dataset
#' #yft.predC <- predict(yft.pr, type = "class")   # predicted 
#' #classification
#'          
#' @importFrom "stats" "na.pass" ".checkMFClasses" "model.frame" "na.omit" "var"                                         
#' @export
predict <- function(object, ...){
      UseMethod("predict")
}

#' @rdname predict
#' @export
predict.dpart <- function(object, newdata = list(), type = c("prob", "class"), 
                        node.cols = NULL,
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
        newdata <- model.frame(Terms, newdata, na.action = na.pass, 
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
        
      #  if(dev.cur() == 1) windows(8,8, record = TRUE)
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
     #   if(dev.cur() == 1) windows(8,8, record = TRUE)
        par(mfrow=c(2,2))
        for(i in names(table(where))){
 
          loss <- object$frame$dev[as.integer(i)]/object$frame$wt[as.integer(i)]
          explore(object = object, pred=pred, loss = loss, pred.where=where,
                  node=i, showtitle = TRUE, cex = cex, cols = node.cols)
        }
          
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




