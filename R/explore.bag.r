explore.bag <- function(object, node, cols = NULL, showtitle = FALSE, axis.side = 2, cex = 1.0, ylim){
  # test
   if(!(axis.side == 2 | axis.side ==4))
      stop("Incorrect axis.side specified. Only 2 or 4 accepted.")
   nc <- ncol(object$m)
   if(missing(ylim))
      ylimit <- c(-0.05, 1.05)
   else
      ylimit <- ylim
   bID <- match(node, row.names(object$m))

   bars <- barplot(object$m[bID,], plot = FALSE)
   x <- as.vector(bars)

   
   nodevals <- data.frame(m = object$m[bID,], lci95 = object$lci95[bID,], 
                          uci95 = object$uci95[bID,], prey = names(data.frame(object$m)))
   tmp <- as.matrix(nodevals[,2:3])
   id <- apply(tmp, 1, function(x){ all(x == 0)})
   nodevals[id,1:3] <- NA
   
   preyO <- 1:length(nodevals$prey)
   
   
   p <- ggplot(nodevals, mapping = aes(x = reorder(prey, preyO), y = lci95)) + 
     geom_segment(stat = "identity", aes(xend = prey, yend = uci95, 
                                         colour = reorder(cols, preyO)), 
                  lineend = "butt", size = 1.5, 
                  arrow = arrow(ends = "both", angle = 90, length = unit(0.1, "cm"),
                                type = "closed")) + 
     scale_color_manual(values = as.vector(cols), labels = names(cols), 
                        name = "Prey")  +
     ylim(ylim) + xlab("") + ylab("Bootstrapped Proportion") + 
     ggtitle(paste("Node", node)) + theme_bw() +
     theme(plot.title = element_text(hjust = 0.5, size = 16),
           plot.margin = unit(c(1,1,1.5,1.2), "cm")) 
   
   print(p)

   data.frame(m = object$m[bID,], v = object$v[bID,], lci95 = object$lci95[bID,],
       uci95 = object$uci95[bID,])
}




