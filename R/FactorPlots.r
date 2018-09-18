#' FactorPlots
#' 
#' @description Internal function used by \code{plot.diet}
#' 
#' @param x x
#' @param XvarF XvarF
#' @param i i
#' @param PredSpID PredSpID
#' 
#' @importFrom "ggplot2" "ggplot" "ggplot_gtable" "aes_string" "geom_point" "geom_bar"
#' @importFrom "ggplot2" "geom_histogram"
FactorPlots <- function(x, XvarF, i, PredSpID){
  
  levels(XvarF[,i])[levels(XvarF[,i]) == ""] <- NA
  x.tab <- tapply(XvarF[,i], x[,PredSpID], function(x) table(x))
  n <- unlist(lapply(x.tab, length))
  x.tab.df <- data.frame(cbind(unlist(x.tab), rep(levels(XvarF[,i]), length(unique(x[,PredSpID]))),
                               rep(names(n), n)))
  names(x.tab.df) <- c("y", "X", "PredSpp")
  x.tab.df$y <- as.numeric(as.vector(x.tab.df$y))
  row.names(x.tab.df) <- NULL
  tab <- tapply(x.tab.df$y, x.tab.df$X, function(x) x/sum(x))
  x.bp <- data.frame(cbind(unlist(tab), rep(names(tab), lapply(tab, length)),
                           rep(levels(x[,PredSpID]), length(tab))))
  names(x.bp) <- c("y", "X", "Predator")
  x.bp$y <- as.numeric(as.vector(x.bp$y))  
  expl9 <- ggplot(data = x.bp, aes_string("X", col = "PredSpp")) + geom_bar() +
    xlab(names(XvarF)[i])
  expl9
}