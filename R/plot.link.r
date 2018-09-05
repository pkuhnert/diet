plot.link <- function(x, ...){
  
  
  if (!inherits(x, "link"))
    stop("Not a link object. \n")
  
  
 
 
  def.par <- par()$mar

  # Plot 1
  layout(matrix(c(1, 2), byrow = TRUE, ncol = 2), heights = c(5,5), widths = c(1,4), respect = TRUE)
  par(mar = c(5, 4, 5, 1))
  image.cols <- c("white", rev(heat.colors(19)))
  z <- matrix(seq(0, 1, length = 20), byrow = TRUE, ncol = 1)
  image(1, 1:20, t(z), col = image.cols, axes = FALSE, xlab = "", ylab = "")
  axis(side = 2, labels = round(seq(0, 1, length = 10), 2), at = seq(1, 20, length = 10), las = 2)
  box()
  

  # Plot 2 
  par(mar = c(9,5,3,2)+0.1)
  image(1:ncol(x$m), 1:nrow(x$m), t(x$m), zlim = c(0,1), 
        col = image.cols, axes = FALSE, xlab = "", ylab = "")
  box()
  mtext(side = 2, at = 1:nrow(x$m), row.names(x$m), line = 1, adj = 1, cex = 0.6, las = 2)
  mtext(side = 2, at = nrow(x$m)/2, "Terminal Nodes", line = 3)
  mtext(side = 1, at = 1:ncol(x$m), names(data.frame(x$m)), line = 0.5, adj = 1, cex = 0.8,
       las = 2)
  par(fig = c(0,1,0,1), mar = c(5,4,4,2)+0.1)
  mtext(side = 3, text = "Predicted Prey Composition at Terminal Nodes", line = 2)
  
  invisible(par(def.par))
  
  
  
}