writepn.csv <- function(x){
  
  if (!inherits(x, "diet"))
    stop("Not an object of class diet")
  
  lvls <- length(levels(x$Group))
  dat <- data.frame(Group = sort(levels(x$Group)), PreyTaxonSort = rep("", lvls), 
                    PreyTaxBroad = rep("", lvls), Sort = rep("", lvls))
  write.csv(dat, "PreyTaxonSort_Template.csv", row.names = FALSE)
  cat("Prey groupings written to file: PreyTaxonSort.csv, Please edit file to assign prey
      groupings. \n")
  
  
}