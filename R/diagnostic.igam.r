diagnostic.igam <- function(obj){
  
  par(mfrow=c(2,2))
  # resids vs fitted
  plot(fitted(obj), obj$residuals, pch = 16, xlab = "Fitted Values", ylab = "Residuals",
       main = "Residuals vs Fitted")
  abline(h=0, lty = 3, col = "red")
  # qqnorm 
  qqnorm(obj$residuals, main = "Q-Q Plot of Residuals")
  qqline(obj$residuals)
  # histogram of residuals
  hist(obj$residuals, main = "Histogram of Residuals")
  # observed vs predicted
  xyrange <- range(fitted(obj), obj$y)
  plot(fitted(obj), obj$y, xlab = "Fitted", ylab = "Response", xlim = xyrange, ylim = xyrange,
       main = "Observed vs Predicted", pch = 16)
  abline(0,1, col = "blue")
  par(mfrow=c(1,1))
 
  
  
  
}