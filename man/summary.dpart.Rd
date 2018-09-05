\name{summary.dpart}
\alias{summary.dpart}
\title{
Summarise a Fitted Dpart Object
}

\description{
Outputs a detailed listing of the classification tree to file: tree_output.txt 
}
\usage{
\method{summary}{dpart}(object, cp = 0, digits = getOption("digits"), 
      file = "tree_output.txt", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
Describe \code{x} here
}
  \item{cp}{
trim nodes with a complexity of less than cp from the listing.
}
  \item{digits}{
Number of significant digits to be used in the result.
}
  \item{file}{
write the output to a given file name. (Full listings of a tree are often quite long).
Default file name is \code{tree_output.txt}.
}
  \item{\dots}{
additional arguments affecting the summary produced
}

}
\details{
This function is a method for the generic function summary for class "dpart". 
It can be invoked by calling summary for an object of the appropriate class, 
or directly by calling \code{summary.dpart} regardless of the class of the object. 
By default the summary is written to file.
}


\seealso{
\code{\link{summary}}, \code{\link{dpart}}, \code{\link{print.dpart}} 
}
\examples{
# Load data
data(yftdiet)  
# Load the prey taxa data
data(PreyTaxonSort)

# Assigning prey colours for default palette
val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = TRUE)
node.colsY <- val$cols
dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes

# Fitting the classification tree
yft.dp <- dpart(Group ~ Lat + Lon + Year + Quarter + SST  + Length, 
                   data = dietPP, weights = W, minsplit = 10,
                   cp = 0.001)
plot(yft.dp, node.cols = node.colsY)
summary(yft.dp)
print(yft.dp, setID = "TripSetNo")
}
