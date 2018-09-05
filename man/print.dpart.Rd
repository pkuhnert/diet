\name{print.dpart}
\alias{print.dpart}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print a Dpart Pbkect
}

\description{
This function prints an \code{dpart} object. It is a method for the generic function
\code{print} of class \code{dpart}.
}
\usage{
\method{print}{dpart}(x, setID = NULL, digits=getOption("digits"), 
   file = "diet_tree_summary.csv", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
fitted model object of class \code{dpart}. 
}
  \item{setID}{
the set identification number used to summarise the output. Defaults to NULL.
}

  \item{digits}{
the number of digits of numbers to print
}

  \item{file}{
file name where to output summary results. (Output too big to print to screen).
}

  \item{\dots}{
arguments to be passed to or from other methods
}
}
\details{
This function is a method for the generic function \code{print} for class "dpart".
It can be invoked by calling print for an object of the appropriate class, or
directly by calling print.dpart.

The loss is computed as the deviance at a node divided by the number of predators
appearing at a node and ranges between 0 and 1. The loss is a measure of the 
diversity of prey eaten by predators at the node, where values near 0 indicate
low diversity and values near 1 indicate high diversity.
}
\value{
Summary of the tree is provided as a csv file with the following column headings:
\item{node}{node number}
\item{nobs}{number of observations in each node}
\item{nsets}{number of sets in each node if a setID was provided otherwise, NA}
\item{npredator}{number of predators in each node}
\item{nprey}{number of prey in each node}
\item{dev}{node deviance}
\item{loss}{node loss}
\item{split}{the split that derived each node}
\item{pclass}{the predicted classification for each node}
}



\seealso{
\code{\link{print}}, \code{\link{diversity}}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
