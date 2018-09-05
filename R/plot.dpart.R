plot.dpart <- function(x, node.cols = NULL, keep.margins = FALSE, ...){


  

  # Plot 2: tree
  val <- prp(x, extra = 8, box.palette = as.list(node.cols), 
      shadow.col = "gray", split.prefix = "is ", 
      split.suffix = "?", nn = TRUE, split.cex = 1.5, cex = 0.6, yesno = TRUE,
      type = 2, branch.lwd = 3, legend.x = NULL, legend.y = NULL)
  
 


  invisible(val)
}




