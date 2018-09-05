dpart <-
function(formula, data, weights, subset, na.action = na.dpart,
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
  