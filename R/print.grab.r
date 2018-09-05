print.grab <- function (x, ...) 
{
  
  object <- x
    cat("Node \t", "N \t", "No. of sets \t", "No. of predators \t", 
        "No. Prey \t", "Deviance \t", "Expected Loss \t Predicted Class\n")
    cat("------------------------------------------------------------------------------------------------------------------------------------\n")
    for (i in 1:nrow(object$nodeS))
        with(object$nodeS[i, ], cat(node,
        "\t", nobs, "\t", nsets, "\t\t", npredators, "\t\t\t", 
        nprey, "\t\t", round(dev, 2), "\t\t", round(loss, 4), "\t", as.character(pclass),
        "\n"))

     # Printing to file
     write.csv(object$nodeS, "grab_node_res.csv", row.names = FALSE)
     cat("Written grabbed node summary table out to file: grab_node_res.csv\n")

     # Plotting bootstrap predictions
     if(any(names(object) == "nodestats")){
       nodestats <- object$nodestats
       options(warn = -1)
       par.old <- par()$mar
       par(mar = c(7, 4, 4, 2) + 0.1)
       plot(1:nrow(nodestats), nodestats$m, type = "n", axes = FALSE,
          ylim = c(0,1), xlab = "", ylab = "Bootstrapped Proportion")
       points(1:nrow(nodestats), nodestats$m, pch = 16)
       arrows(1:nrow(nodestats), nodestats$lci95, 1:nrow(nodestats), nodestats$uci95,
         angle = 90, code = 3, length = 0.1)
       box()
       axis(side = 2)
       mtext(side = 1, text = row.names(nodestats), adj = 1, line = 1,
          at = 1:nrow(nodestats), las = 2)

       title(main = paste("Node", object$nodeS$node))
       par(mar = par.old)
       options(warn = 0)
      }
     invisible()
}
