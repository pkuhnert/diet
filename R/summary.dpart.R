#' Summarise a Fitted dpart Object
#' 
#' @description Outputs a detailed listing of the classification tree to file: tree_output.txt 
#' 
#' @param object dpart object
#' @param cp trim nodes with a complexity of less than cp from the listing.
#' @param digits Number of significant digits to be used in the result.
#' @param file write the output to a given file name. (Full listings of a tree are often quite long).
#' Default file name is \code{tree_output.txt}.
#' @param ... additional arguments affecting the summary produced
#' 
#' 
#' @details   This function is a method for the generic function summary for class \code{dpart}. 
#' It can be invoked by calling summary for an object of the appropriate class, 
#' or directly by calling \code{summary.dpart} regardless of the class of the object. 
#' By default the summary is written to file.
#' 
#' @seealso   \code{\link{summary}}, \code{\link{dpart}}, \code{\link{print.dpart}} 
#' 
#' @examples 
#' # Load data
#' #data(yftdiet)  
#' 
#' # Load the prey taxa data
#' #data(PreyTaxonSort)
#' 
#' # Assigning prey colours for default palette
#' #val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = TRUE)
#' #node.colsY <- val$cols
#' #dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes
#' 
#' # Fitting the classification tree
#' #yft.dp <- dpart(Group ~ Lat + Lon + Year + Quarter + SST  + Length, 
#'#                   data = dietPP, weights = W, minsplit = 10,
#' #                                    cp = 0.001)
#'                                     
#' #plot(yft.dp, node.cols = node.colsY)
#' #summary(yft.dp)
#' #print(yft.dp, setID = "TripSetNo")
#' @export
summary <- function(object, ...) 
  UseMethod("summary")

#' @rdname summary
summary.dpart <- function(object, cp = 0, digits = getOption("digits"), file = "tree_output.txt", ...){


 summary.rpart(object, cp = cp, digits = digits, file = file, ...)
 wd <- getwd()
 cat(paste("Summary information written to: ", wd, "/", file, "\n", sep = ""))

}


#' @rdname summary
#' @import "stats"
summary.rpart <- function(object, cp=0, digits=getOption("digits"), file,  ...)
  {
    if(!inherits(object, "rpart")) stop("Not legitimate rpart object")
    
    # If this is an older-style rpart object, convert it
    #  either way, rename it to "x" to save typing (commented this out as compiler didn't like it)
    # if (!is.null(object$frame$splits)) x <- rpconvert(object)
    # else  x <- object
    
    x <- object
    
    if (!missing(file)) {
      sink(file)
      on.exit(sink())
    }
    
    if(!is.null(x$call)) {
      cat("Call:\n")
      dput(x$call, control=NULL)
    }
    
    omit <- x$na.action
    n <- x$frame$n
    if (length(omit))
      cat("  n=", n[1L], " (", naprint(omit), ")\n\n", sep="")
    else cat("  n=", n[1L], "\n\n")
    
    print(x$cptable, digits=digits)
    ff <- x$frame
    ylevel <- attr(x,'ylevels')
    id <- as.integer(row.names(ff))
    parent.id <- ifelse(id==1,1, floor(id/2))
    parent.cp <- ff$complexity[match(parent.id, id)]
    rows <- (1L:length(id))[parent.cp > cp]
    if (length(rows)>0) rows <- rows[order(id[rows])]
    else rows <- 1L
    is.leaf <- (ff$var=='<leaf>')
    index <- cumsum(c(1, ff$ncompete + ff$nsurrogate + 1*(!is.leaf)))
    
    if(!all(is.leaf)) {  #skip these lines for a "no splits" tree
      sname <- dimnames(x$splits)[[1L]]
      cuts <- vector(mode='character', length=nrow(x$splits))
      temp <- x$splits[ ,2L]
      for (i in 1L:length(cuts)) {
        if(temp[i] !=0){ 
          if (temp[i] == -1L)
            cuts[i] <-paste("<", format(signif(x$splits[i,4], digits=digits)))
          else if (temp[i] ==1L)
            cuts[i] <-paste("<", format(signif(x$splits[i,4], digits=digits)))
          else cuts[i]<- paste("splits as ",
                               paste(c("L", "-", "R")[x$csplit[x$splits[i,4], 1L:temp[i]]],
                                     collapse='', sep=''), collapse='')
        }
      }
      # S-PLUS 4.0 can't handle null vectors here
      if(any(temp<2)) cuts[temp<2 ] <- format(cuts[temp<2],justify="left")
      cuts <- paste(cuts, ifelse(temp >=2, ",",
                                 ifelse(temp==1, " to the right,", " to the left, ")),
                    sep = '')
    }
    
    if (is.null(ff$yval2))
      tprint <- x$functions$summary(ff$yval[rows], ff$dev[rows],
                                    ff$wt[rows], ylevel, digits)
    else
      tprint <- x$functions$summary(ff$yval2[rows,,drop=TRUE], ff$dev[rows],
                                    ff$wt[rows], ylevel, digits)
    
    for (ii in 1L:length(rows)) {
      i <- rows[ii]
      nn <- ff$n[i]
      cat("\nNode number ", id[i], ": ", nn, " observations", sep='')
      if (ff$complexity[i] < cp || is.leaf[i]) cat("\n")
      else cat(",    complexity param=",
               format(signif(ff$complexity[i], digits)), "\n", sep="")
      
      cat(tprint[ii], "\n")
      if (ff$complexity[i] > cp && !is.leaf[i] ){
        sons <- 2*id[i] + c(0,1)
        sons.n <- ff$n[match(sons, id)]
        cat("  left son=", sons[1L], " (", sons.n[1L], " obs)",
            " right son=", sons[2L], " (", sons.n[2L], " obs)", sep='')
        j <- nn - (sons.n[1L] + sons.n[2L])
        if (j>1L) cat(", ", j, " observations remain\n", sep='')
        else if (j==1L) cat(", 1 observation remains\n")
        else     cat("\n")
        cat("  Primary splits:\n")
        j <- seq(index[i], length.out=1L+ff$ncompete[i])
        if (all(nchar(cuts[j], "w") < 25))
          temp <- format(cuts[j], justify="left")
        else  temp <- cuts[j]
        cat(paste("      ", format(sname[j], justify="left"), " ", temp,
                  " improve=", format(signif(x$splits[j,3], digits)),
                  ", (", nn - x$splits[j,1L], " missing)", sep=''),
            sep="\n")
        if (ff$nsurrogate[i] >0L) {
          cat("  Surrogate splits:\n")
          j <- seq(1 +index[i] + ff$ncompete[i], length.out=ff$nsurrogate[i])
          agree <- x$splits[j,3L]
          if (all(nchar(cuts[j], "w") < 25))
            temp <- format(cuts[j], justify="left")
          else  temp <- cuts[j]
          if (ncol(x$splits)==5) {
            adj   <- x$splits[j,5L]
            cat(paste("      ", format(sname[j], justify="left"), " ",
                      temp,
                      " agree=", format(round(agree, 3)),
                      ", adj=" , format(round(adj, 3)),
                      ", (", x$splits[j,1L], " split)", sep=''),
                sep="\n")
          }
          else {                  #an older style rpart object -- no adj value present
            cat(paste("      ", format(sname[j], justify="left"), " ",
                      temp,
                      " agree=", format(round(agree, 3)),
                      ", (", x$splits[j,1L], " split)", sep=''),
                sep="\n")
          }
        }
      }
    }
    cat("\n")
    invisible(x)
  }


