\name{diversity.dpart}
\alias{diversity}
\alias{diversity.dpart}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diversity Plots
}

\description{
This function produces maps of diversity based on a fitted tree object of class 
\code{dpart}.
}
\usage{
\method{diversity}{dpart}(object, newdata = NULL, LatID, LonID,
     mapxlim = NULL, mapylim = NULL, cex.axis = 0.75,
     mapcol = "gold3", plot = TRUE, database = "world")
}
\arguments{
  \item{object}{
fitted tree object of class \code{dpart}.
}
  \item{newdata}{
optional new dataset to base the diversity estimates on. By default, the dataset
used to create the tree is used.
}
  \item{LatID}{
column labelling the latitude co-ordinate
}
  \item{LonID}{
column labelling the longitude co-ordinate
}
  \item{mapxlim}{
x-range for the map. Defaults to the range of the longitude co-ordinate.
}
  \item{mapylim}{
y-range for the map. Defaults to the range of the latitude co-ordinate.

}
  \item{cex.axis}{
size of points on both axes. Defaults to 0.75.
}
  \item{mapcol}{
colour of land masses shown on the map. Defaults to "gold3"
}
  \item{plot}{
logical. If true, a map is produced.
}
  \item{database}{
one of "world" or "world2" specifying the map co-ordinate system to use. If no
database is given, this will default to "world".
}
  
}
\details{
Depending on the extent of the data collected and co-ordinate system required, 
"world2" may produce a better map than "world". 
}
\value{
Two maps showing the diversity of predators diet spanning the sampled locations.
The first map shows the diversity for each sampled point represented by a colour
matched to a 0-1 scale. The second plot shows the interpolated diversity from
fitting a GAM (see \code{\link{mgcv}}) to the diversity score with smooth terms
for latitude and longitude.

A data frame listing the latitude, longitude, diversity, diversity bin and associated
colours.
}
\references{
Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
and Stable Isotope Data, Journal of Statistical Software, In Prep. 

}
\author{
Petra Kuhnert and Leanne Duffy
}



\seealso{
\code{\link{mgcv}}, \code{\link{mappoints.data}}, \code{\link{map}}
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
diversity(object = yft.pr, LatID = "Lat", LonID = "Lon", too.far = 0.05,
          projection = "+proj=longlat +datum=WGS84", maxpixels = 1e5,
          col = "black", fill = "lightgrey")
}
