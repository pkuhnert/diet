#' Diet Partitioning
#' 
#' @description Analyses diet data using a classification tree analysis with case weights
#' 
#' @usage dpart(formula, data, weights, subset, na.action = na.dpart,
#' model = FALSE, x = FALSE, y = TRUE, parms, control, cost, ...)
#' 
#' @param formula  a \link{formula}, with a response but no interaction terms as 
#' for the \code{\link{rpart}} function
#' @param data an optional data frame in which to interpret the variables named 
#' in the formula
#' @param weights case weights
#' @param subset optional expression saying that only a subset of the rows of the 
#' data should be used in the fit.
#' @param na.action The default action deletes all observations for which \code{y} 
#' is missing, but keeps those in which one or more predictors are missing.
#' @param model if logical: keep a copy of the model frame in the result? If the input value for
#' model is a model frame (likely from an earlier call to the rpart function), then 
#' this frame is used rather than constructing new data.
#' @param x keep a copy of the x matrix in the result.
#' @param y keep a copy of the dependent variable in the result. If missing 
#' and \code{model} is supplied this defaults to \code{FALSE}.
#' @param parms optional parameters for the splitting function. For classification splitting, 
#' the list can contain any of: the vector of prior probabilities (component prior), 
#' the loss matrix (component loss) or the splitting index (component split). The priors 
#' must be positive and sum to 1. The loss matrix must have zeros on the diagonal and 
#' positive off-diagonal elements. The splitting index can be gini or information. The 
#' default priors are proportional to the data counts, the losses default to 1, and the 
#' split defaults to \code{gini}.
#' @param control options that control details of the \code{rpart} algorithm.
#' @param cost a vector of non-negative costs, one for each variable in the model. Defaults 
#' to one for all variables. These are scalings to be applied when considering splits, 
#' so the improvement on splitting on a variable is divided by its cost in deciding which split to choose.
#' @param \dots arguments to \code{rpart.control} may also be specified in the call to 
#' \code{rpart}. They are checked against the list of valid arguments.
#'
#' @details Analyses diet data using a univariate tree analysis with case weights
#' 
#' @value an object of class \code{dpart}, a superset of class \code{rpart}.
#' 
#' @references 
#' Breiman, Friedman, Olshen, and Stone. (1984) Classification and Regression Trees. Wadsworth.
#' 
#' Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet and Stable
#' Isotope Data, Journal of Statistical Software, In Prep. 
#' 
#' @seealso \link{rpart}
#' 
#' @keywords classification, diet, tree
#' 
#' @examples 
# Load data
data(yftdiet)  
# Load the prey taxa data
data(PreyTaxonSort)
  
# Assigning prey colours for default palette
val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = TRUE)
node.colsY <- val$cols
dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes
  
# Fitting the classification tree
yft.dp <- dpart(Group ~ Lat + Lon + Year + Quarter + SST  + Length, 
                  data = dietPP, weights = W, minsplit = 10,
                  cp = 0.001)
plot(yft.dp, node.cols = node.colsY)
summary(yft.dp)
print(yft.dp, setID = "TripSetNo")



dpart <- function(formula, data, weights, subset, na.action = na.dpart,
   model = FALSE, x = FALSE, y = TRUE, parms, control, cost,
   ...){


    m <- match.call(expand.dots = FALSE)
    m$model <- m$method <- m$control <- NULL
    m$x <- m$y <- m$parms <- m$... <- NULL
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
    
    if(!is.factor(Y))
      stop("Response needs to be a factor e.g. benthic classification for this type of analysis.\n")
    else{
      method <- "class"
      if(missing(subset) & missing(parms))
        val <- rpart(formula, data, weights = W, na.action = na.action,
                     method = method, model = model, x = x, y = y, control = control, cost = cost,
                     ...)
      else if(missing(subset)){
          val <- rpart(formula, data, weights = W, na.action = na.action,
                       method = method, model = model, x = x, y = y, parms = parms, control = control,
                       cost = cost,  ...)
      }
      else if(missing(parms)){ 
        val <- rpart(formula, data, weights = W, subset, na.action = na.action,
                     method = method, model = model, x = x, y = y, control = control, cost = cost,
                     ...)
      }
      else
        val <- rpart(formula, data, weights = W, subset, na.action = na.action,
                     method = method, model = model, x = x, y = y, parms, control = control,
                     cost = cost,  ...)
    }
    
      
    
    
    val$data <- data
    oldClass(val) <- c("dpart", class(val))
    val
    
  }
  