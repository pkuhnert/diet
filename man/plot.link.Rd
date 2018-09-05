\name{plot.link}
\alias{plot.link}

\title{
Plots Prey Composition
}
\description{
Produces an image plot of the prey composition at each terminal node of the tree using
the bagged estimates
}
\usage{
\method{plot}{link}(x, ...)
}
\arguments{
  \item{x}{
object of class \code{link}.
}
  \item{\dots}{
arguments to be passed to or from other methods
}
}


\references{
Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
and Stable Isotope Data, Journal of Statistical Software, In Prep. 
}
\author{
Petra Kuhnert and Leanne Duffy
}

\seealso{
\code{\link{image}}
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
yft.pr <- prune(yft.dp, se = 1)
plot(yft.pr, node.cols = node.colsY)


# Bagging
# Bagging with NO spatial bootstrapping
yft.bag <- bagging(Group ~ Lat + Lon + Year + Quarter + SST  + Length,
                    data = dietPP, weights = W, minsplit = 50,
                    cp = 0.001, nBaggs = 500, predID = "TripSetPredNo")

# Link the predictions back to terminal nodes of a tree and plot 
ex.bag.l <- link(x = yft.bag, object = yft.pr, LatID = "Lat", LonID = "Lon",
                 mapxlim = c(-125, -75), mapylim = c(0, 30), plot = TRUE, oob = TRUE, 
                 mfrow = c(2,2))
 
# Plot prey predictions
plot(ex.bag.l)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
