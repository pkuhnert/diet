explore <- function(object, pred, pred.where, loss = NULL, node, cols = NULL,
   showtitle = FALSE, labels = TRUE, cex = 1.0, ylim){

  pred.node <- pred[pred.where == paste(node),]
  if(length(pred.node) == 0)
    print(cat(paste("No data at node", node, "\n")))
  else if(is.null(nrow(pred.node)))
    pred.node.m <- pred.node
  else
    pred.node.m <- apply(pred.node, 2, mean)

  if(showtitle)
     par(mar = c(6,4,4,2)+0.1)
  n <- length(pred.node.m)

  bars <- barplot(pred.node.m, plot = FALSE)
  x <- as.vector(bars)
  if(labels)
     nms <- names(pred.node.m)
  else
     nms <- NULL
  names(pred.node.m) <- NULL
  if(missing(cols))
     cols <- "grey"
  
  if(missing(ylim))
    ylim <- c(-0.05,1.05)


  preyO <- 1:length(pred.node)
  p <- ggplot(mapping = aes(x = reorder(names(pred.node),preyO), 
                       y = pred.node.m)) + 
    geom_bar(stat = "identity", aes(fill = reorder(cols, preyO))) + 
    scale_fill_manual(values = as.vector(cols), labels = names(cols), name = "Prey")  +
    ylim(ylim) + xlab("") + ylab("Proportion") + 
    ggtitle(paste("Node ", row.names(object$frame)[as.integer(node)], "\n",
            "Diet Composition (D=", round(loss, 3), ")", sep = "")) + theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 16),
          plot.margin = unit(c(1,1,1.5,1.2), "cm")) 
  
  print(p)
 

}