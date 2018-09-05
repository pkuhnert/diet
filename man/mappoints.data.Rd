\name{mappoints}
\alias{mappoints}
\alias{mappoints.data}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Produce a Map of a Region with Data overlayed
}

\description{
Utilises the \code{map} function to produce a map of the study area based on
locations given.
}
\usage{
\method{mappoints}{data}(XY.data, xlim, ylim, mar, add = FALSE, col = "gold3", 
ptcol = "grey", database = 'world', cex.axis = 0.75, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XY.data}{
matrix of longitude and latitude co-ordinates
}
  \item{xlim}{
map x limits
}
  \item{ylim}{
map y limits
}
  \item{mar}{
optional margins used for producing the map
}
  \item{add}{
logical. Whether to add the map to an existing plot.
}
  \item{col}{
one of \code{colors()} that indicates the colour for the land mass plotted on
the map
}
  \item{ptcol}{
colour of the points plotted on the map. Default = "grey"
}
  \item{database}{
either "world" or "world2" depending on the area being plotted. Defaults to "world"
}
  \item{cex.axis}{
axis label size. Default = 0.75.
}
  \item{\dots}{
Extra arguments passed to the \code{\link{map}} function.
}
}
\details{
Calls on the \code{map} function for plotting.
}

\references{
Richard A. Becker, and Allan R. Wilks, "Maps in S", AT&T Bell Laboratories 
Statistics Research Report [93.2], 1993. 
http://public.research.att.com/areas/stat/doc/93.2.ps

Richard A. Becker, and Allan R. Wilks, "Constructing a Geographical Database", 
AT&T Bell Laboratories Statistics Research Report [95.2], 1995. 
http://public.research.att.com/areas/stat/doc/95.2.ps

}


\seealso{
\code{\link{map}},
}
\examples{
# Load dataset
data(yftdiet)

# Produce map of the data
mappoints(XY.data = yftdiet[,c("Lon", "Lat")])
title(main = "YFT Data Distribution")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
