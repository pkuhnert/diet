summary.dpart <-
function(object, cp = 0, digits = getOption("digits"), file = "tree_output.txt", ...){


 summary.rpart(object, cp = cp, digits = digits, file = file, ...)
 wd <- getwd()
 cat(paste("Summary information written to: ", wd, "/", file, "\n", sep = ""))

}

