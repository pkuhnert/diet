\name{prune.dpart}
\alias{prune.dpart}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cost-complexity pruning of a Dpart object
}

\description{
This function prunes a fitted \code{dpart} object based on either the standard
error or number of splits required.
}
\usage{
\method{prune}{dpart}(tree, se, nsplits, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
a legitimate tree object of class \code{dpart}.
}
 \item{se}{
numeric. A standard error used to prune. 
}
 \item{nsplits}{
numeric. Number of splits to prune to.
}
  \item{\dots}{
arguments to be passed to or from other methods
}
}

\value{
A new \code{dpart} object that is pruned based on the number of standard errors
or number of splits.
}
\references{
Breiman, L., Friedman, J.H., Olshen, R.A. and Stone, C.J. (1984) Classification 
and Regression Trees. Wadsworth International.
}
\author{
Petra Kuhnert
}



\seealso{
\code{\link{dpart}}
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
