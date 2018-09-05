igam <- function(formula, data, control = list(), G = NULL, ...){
    
  
   
   gp <- interpret.gam(formula)
   m <- match.call(expand.dots = FALSE)
   m$formula <- gp$fake.formula
   Y <- gp$response

   
   m$model <- m$method <- m$control <- NULL
   m$x <- m$y <- m$parms <- m$... <- NULL
#   m$na.action <- na.action
   m$family <- m$control <- m$scale <- m$knots <- m$sp <- m$min.sp <- m$H <- m$select <- m$gamma <- m$method <- m$fit <- m$paraPen <- m$G <- m$optimizer <- m$in.out <- m$... <- NULL
   m$drop.unused.levels <- TRUE
   m[[1]] <- as.name("model.frame")
   m <- eval(m, parent.frame())

    val <- gam(formula=formula, family = gaussian(), data = data,  ...)


    val$data <- data
    oldClass(val) <- c("igam", class(val))

 val
}