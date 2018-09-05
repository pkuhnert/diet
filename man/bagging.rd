\name{bagging}
\alias{bagging}

\title{
Bagging 
}

\description{
Creates bagged tree estimates from diet data
}
\usage{
bagging(formula, data, weights, subset, na.action = na.dpart,
      model = FALSE, x = FALSE, y = TRUE, parms, control, cost, nBaggs,
      spatial = list(fit = FALSE, sizeofgrid = 5, nsub = NULL, ID = NULL, LonID = "Longitude", 
                     LatID = "Latitude"), Plot = FALSE,
       predID, numCores = 1, ...)

}
\arguments{
  \item{formula}{
a \link{formula}, with a response but no interaction terms as for the \code{\link{rpart}}
 function
}
  \item{data}{
an optional data frame in which to interpret the variables named in the formula

}
  \item{weights}{
case weights

}
  \item{subset}{
optional expression saying that only a subset of the rows of the data should be used in the fit.
}
  \item{na.action}{
The default action deletes all observations for which \code{y} is missing, but keeps
those in which one or more predictors are missing.
}
  \item{model}{
if logical: keep a copy of the model frame in the result? If the input value for
model is a model frame (likely from an earlier call to the rpart function), then
this frame is used rather than constructing new data.
}
  \item{x}{
keep a copy of the x matrix in the result.
}
  \item{y}{
keep a copy of the dependent variable in the result. If missing and \code{model}
is supplied this defaults to \code{FALSE}.
}
  \item{parms}{
optional parameters for the splitting function. 
For classification splitting, the list can contain any of: the
vector of prior probabilities (component prior), the loss matrix (component loss) or
the splitting index (component split). The priors must be positive and sum to 1.
The loss matrix must have zeros on the diagonal and positive off-diagonal elements.
The splitting index can be gini or information. The default priors are proportional
to the data counts, the losses default to 1, and the split defaults to \code{gini}.
}
  \item{control}{
options that control details of the \code{rpart} algorithm.
}
  \item{cost}{
a vector of non-negative costs, one for each variable in the model.
Defaults to one for all variables. These are scalings to be applied when
considering splits, so the improvement on splitting on a variable is divided by
its cost in deciding which split to choose.
}
  \item{nBaggs}{
numeric. Number of bootstrap samples.
}
  \item{spatial}{
A list with the following elements:
         fit  = do spatial bootstrapping
         sizeofgrid =  size of spatial tile to sample from (default is 5)
         nsub = number of sub-samples to take (defaults to no subsampling)
         ID   = ID in which to subsample from (e.g. TripSetPredNo)
                (only required if sub-sampling is required)
}
  \item{Plot}{
plotting the spatial grid with samples (default: no plotting (FALSE))
}   \item{predID}{
predator ID
}
 \item{numCores}{
Number of cores to push the bagging on to. Only available under Unix (default: 1)
}
  \item{\dots}{
arguments to be passed to or from other methods
}
}
\details{
Users will need to determine whether spatial bootstrapping is required. They can
use the \code{\link{resid}} function to examine the residuals from the fit of the
model to determine whether this is required.
}
\value{
A list with the following elements:
\item{baggs}{tree objects for each \code{B} trees produced.}
\item{oob}{numeric vector indicating the samples left as out of bag (oob) samples.}
\item{pred.oob}{predicted prey composition for each set of out of bag samples.}
\item{pred}{all predicted prey compositions for each bootstrap sample.}
\item{resid}{data frame of residuals from the fitted tree for each bootstrap sample.}
\item{data}{bootstrap sample dataset}
}
\references{
Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
and Stable Isotope Data, Journal of Statistical Software, In Prep. 

Kuhnert PM, Kinsey-Henderson A, Bartley R, Herr A (2010)
Incorporating uncertainty in gully erosion calculations using
the random forests modelling approach. Environmetrics
21:493-509. doi:10.1002/env.999

Breiman L (1996) Bagging predictors. Mach Learn 24:123-140. doi:
10.1023/A:1018054314350

Breiman L (1998) Arcing classifiers (with discussion). Ann Stat
26:801-824. doi:10.2307/120055

Breiman L (2001) Random forests. Mach Learn 45:5-32. doi:
10.1023/A:1010933404324

}
\author{
Petra Kuhnert and Leanne Duffy
}



\seealso{
\code{\link{dpart}}; \code{\link{link.bag}} 
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

# Bagging
# Bagging with NO spatial bootstrapping
yft.bag <- bagging(Group ~ Lat + Lon + Year + Quarter + SST  + Length,
                    data = dietPP, weights = W, minsplit = 50,
                    cp = 0.001, nBaggs = 500, predID = "TripSetPredNo")


}
\keyword{ ~bootstrap }
\keyword{ ~sampling }% __ONLY ONE__ keyword per line
\keyword{ ~uncertainty }% __ONLY ONE__ keyword per line
