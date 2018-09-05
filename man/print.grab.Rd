\name{print.grab}
\alias{print.grab}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prints Output from a Node Selected from a Tree
}

\description{
Produces a summary of the node selected using the \code{grab} function.
}
\usage{
\method{print}{grab}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object of class \code{grab}.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}

\value{
Summarises the node selected into the following:
Node ID, Number of observations, Number of sets, Number of predators,
Number of prey, deviance, expected loss and predicted class.

The output is printed to a .csv file
}
\references{
Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
and Stable Isotope Data, Journal of Statistical Software, In Prep. 

}
\author{
Petra Kuhnert
}


\seealso{
\code{\link{print}},
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
                   
# Pruning the tree
yft.pr <- prune(yft.dp, se = 1)                   

# Exploring Nodes: This suite of graphics is interactive and therefore has
# been commented out. When run, the code will ask you to select a node for
# viewing.

# Exploring nodes of the tree - single page
# val <- grab(object = yft.pr, LatID = "Lat", LonID = "Lon", setID = "TripSetNo", 
#            node.cols = node.colsY, cex = 1, mapxlim = c(-125, -75), mapylim = c(0, 30),
#            mapcol = "gold3", pos = "topleft")
# val
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
