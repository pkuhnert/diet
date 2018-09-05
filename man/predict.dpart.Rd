\name{predict.dpart}
\alias{predict}
\alias{predict.dpart}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict from a Dpart Object
}

\description{
Forms predictions from an object of class \code{dpart}. Options are to produce predicted
probabilities or predicted classificaitons.
}
\usage{
\method{predict}{dpart}(object, newdata = list(), type = c("prob", "class"), na.action = na.pass,
   plot = TRUE, pred.type = "obs", predatorID = NULL, cex = 1.0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
object of class \code{dpart}
}
  \item{newdata}{
optional data set to predict on.
}
  \item{type}{
either "prob" or "class" to specify probabilities or classifications.
}
  \item{na.action}{
what to do about missing values. (Default is to omit).
}
  \item{plot}{
logical. Produce a plot of the predictions if TRUE otherwise just write the predictions
out to file.
}
  \item{pred.type}{
either "obs" or "predator" to indicate whether predictions are to be generated for
each observation or each predator, respectively.
}
  \item{predatorID}{
predator identification label. Required when \code{pred.type = "predator"}
}
 \item{cex}{
numeric. Size of plotting characters and labels.
}
  \item{\dots}{
other arguments to pass into the function
}
}

\value{
Predicted probabilities for each observation or each predator, if 
\code{pred.type = "predator"} was specified.
Predicted classifications if \code{type = "class"} was specified.

}
\references{
Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
and Stable Isotope Data, Journal of Statistical Software, In Prep. 

}
\author{
Petra Kuhnert and Leanne Duffy
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

# Predictions
# predict distribution of prey composition for each predator
yft.predator <- predict(yft.pr, type = "prob", pred.type = "predator",
                              predatorID = "TripSetNo")

# predict distribution of prey composition for each observation
yft.pred.obs <- predict(yft.pr, type = "prob")

# predict classification  for each observation in the dataset
yft.predC <- predict(yft.pr, type = "class")   # predicted classification

                   

}
