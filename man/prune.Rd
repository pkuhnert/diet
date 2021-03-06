% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prune.r
\name{prune}
\alias{prune}
\alias{prune.dpart}
\title{Cost-complexity pruning of a dpart object}
\usage{
\method{prune}{dpart}(tree, se, nsplits, ...)

\method{prune}{dpart}(tree, se, nsplits, ...)
}
\arguments{
\item{tree}{a legitimate tree object of class \code{dpart}.}

\item{\dots}{arguments to be passed to or from other methods}

\item{se}{numeric. A standard error used to prune.}

\item{nsplits}{numeric. Number of splits to prune to.}
}
\description{
This function prunes a fitted \code{dpart} object based on either the standard
error or number of splits required.
}
\details{
A new \code{dpart} object that is pruned based on the number of standard errors
or number of splits.
}
\examples{

# Load data
#data(yftdiet)  

# Load the prey taxa data
#data(PreyTaxonSort)

# Assigning prey colours for default palette
#val <- apc(x = yftdiet, preyfile = PreyTaxonSort, check = TRUE)
#node.colsY <- val$cols
#dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes

# Fitting the classification tree
#yft.dp <- dpart(Group ~ Lat + Lon + Year + Quarter + SST  + Length, 
#                   data = dietPP, weights = W, minsplit = 10,
#                                     cp = 0.001)
#yft.pr <- prune(yft.dp, se = 1)
#plot(yft.pr, node.cols = node.colsY)
}
\references{
Breiman, L., Friedman, J.H., Olshen, R.A. and Stone, C.J. (1984) Classification 
and Regression Trees. Wadsworth International.
}
\seealso{
\code{\link{dpart}}
}
