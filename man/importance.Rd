\name{importance.dpart}
\alias{importance}
\alias{importance.dpart}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variable Importance Ranking
}

\description{
Produces a variable importance ranking from a fitted \code{dpart} object.
}
\usage{
\method{importance}{dpart}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object of class \code{dpart}
}

}
\details{
A barplot of the variable importance ranking sorted from highest to lowest is produced 
in addition to the printed output.
}
\value{
Variable importance ranking for the variables used in the classification tree.
}
\references{
Breiman, L., Friedman, J.H., Olshen, R.A. and Stone, C.J. (1984) Classification 
and Regression Trees. Wadsworth International.
}
\author{
Petra Kuhnert
}



\seealso{
\code{\link{dpart}},
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
vi <- importance(yft.pr)
}
