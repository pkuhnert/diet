\name{plot.igam}
\alias{plot.igam}

\title{
Term Plots from GAM
}

\description{
This function is based on the \code{plot.gam} function from the \code{mgcv} package
but with more functionality and applicability for isotope analyses.
}
\usage{
\method{plot}{igam}(x, Xvar, var.cond = NULL, plot.map = FALSE,
                       database = "world", leg.pos = "topleft", 
                       spgrp = NULL, too.far = 0.1, se.plot = TRUE, zlim = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object of class \code{igam} resulting from the GAM fit.
}
  \item{Xvar}{
vector of x-variables to plot. If one x-variable is listed, then a 1D termplot is 
produced. If two x-variables are listed (e.g. c("Lon", "Lat")) then a 2D image plot
is produced. 
}
  \item{var.cond}{
variable name listing the variable to condition on. This will produce multiple plots
according to each unique value in the \code{var.cond} variable.
}

  \item{plot.map}{
logical. Add a map to the set of points. (Default = FALSE)
}
  \item{database}{
map to use for mapping. Either 'world' or 'world2' accepted. Default = 'world'.
See \code{\link{map}} for more details.
}
  \item{leg.pos}{
legend position indicating the position of the legend for 1D plots. See \code{\link{legend}}
for more details of how to position legends. Default: "topleft"
}

  \item{spgrp}{
data frame with two columns named "Group" and Xvar where Xvar is a factor and represents the name of the variable
to be plotted. This dataframe will be used to group the Xvar levels according to the
group assigned.
}

 \item{too.far}{
If greater than 0 then this is used to determine when a location is too far from 
data to be plotted when plotting 2-D smooths. This is useful since smooths tend 
to go wild away from data. The data are scaled into the unit square before deciding 
what to exclude, and too.far is a distance within the unit square. See \code{\link{plot.gam}}
for more details.
}

 \item{se.plot}{
logical. Produce contours of the standard errors (TRUE) otherwise the contours
of the predictions are overlayed.
}

 \item{zlim}{
Z limits
}

  \item{\dots}{
Extra arguments passed to other functions.
}
}

\value{
Linear Terms (Factors): A 1D plot of the predicted response for the \code{Xvar}
specified while holding all other variables constant at their median (Includes the 
constant).

Smooth 1D Terms: A 1D plot of the term contribution in the model for the \code{Xvar}
specified. See \code{\link{mgcv}} for more details.

Smooth 2D Terms: A 2D plot of the predicted response for the \code{Xvars} specified,
while holding all other variables constant at their median.


}
\references{
Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
and Stable Isotope Data, Journal of Statistical Software, In Prep. 

Wood, S.N. (2006) Generalized Additive Models: an introduction with R, CRC
}
\author{
Petra Kuhnert
}


\seealso{
\code{\link{plot.gam}},
}
\examples{
# Load Isotope Dataset
#data(yftISO)

# Fit GAM
#Yd15N.gam <- igam(d15N ~ Year + s(FL_mm) + s(Lon, Lat),
#                     data = yftISO)
#summary(Yd15N.gam)

# Termplots and predictions
# FL_mm
#plot(Yd15N.gam, Xvar = "FL_mm")
# Lat & Lon
#plot(Yd15N.gam, Xvar = c("Lon", "Lat"), database = "world", plot.map = TRUE)
#plot(Yd15N.gam, Xvar = c("Lon", "Lat"), database = "world", too.far = 0.15, 
#  plot.map = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
