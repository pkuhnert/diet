\name{diagnostic}
\alias{diagnostic}
\alias{diagnostic.igam}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostics from Fitted GAM 
}

\description{
Produces a series of 4 diagnostic plots to examine the fit of the GAM model
}
\usage{
\method{diagnostic}{igam}(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
fitted gam object of class \code{igam}
}

}

\value{
Produces 4 plots which can be used to check Normal assumptions:
(1) Residuals versus fitted
(2) Q-Q Normal plot of the Residuals
(3) Histogram of the Residuals
(4) Observed versus Fitted
}
\references{
Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
and Stable Isotope Data, Journal of Statistical Software, In Prep. 
}
\author{
Petra Kuhnert
}



\examples{
# Load Isotope Dataset
#data(yftISO)

# Fit GAM
#Yd15N.gam <- igam(d15N ~ Year + s(FL_mm) + s(Lon, Lat),
#                     data = yftISO)
#summary(Yd15N.gam)

# Check fit of the model
#diagnostic(Yd15N.gam)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
