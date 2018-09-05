\name{apc}
\alias{apc}
\alias{apc.diet}


\title{
Assigns Prey Colours
}

\description{
Assigns prey colours to a matrix of dietmatrix class. 
}
\usage{
\method{apc}{diet}(x, preyfile, palette = NULL, check = FALSE)
}
\arguments{
  \item{x}{
diet object. Must be of class diet.
}
  \item{preyfile}{
prey file name that contains the prey groupings and labels required for labelling
the tree. The preyfile can be output using the \code{writepn.csv} function and then
edited accordingly. See \code{writepn.csv} for more details.
}
\item{palette}{
a vector of colours. This requires the user to set up a colour 
palette where the number of colours equals the number of prey in \code{x}. See
\code{colors()} for a listing of colours available. By default this is set to 
\code{NULL}.
}
\item{check}{
logical. If set to \code{TRUE}, the prey colour palette is plotted.
}
}
\details{
If a colour palette is not provided, then a palette will be automatically 
derived based on the colour palettes available. For customised colour palettes,
users will need to define a \code{palette}. If no palette is supplied and there 
are more than 4 prey groups, the colour palette will default to rainbow.
}
\value{
A list with the following components:
\item{cols}{a vector of node colours for each prey appearing in the diet matrix.}
\item{x}{the updated diet matrix with prey groupings as specified in the 
 PreyTaxonSort.csv file.}

}
\references{
Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
and Stable Isotope Data, Journal of Statistical Software, In Prep. 
}
\author{
Petra Kuhnert and Leanne Duffy
}


\seealso{
\code{\link{colors}}; \code{\link{palette}}; \code{\link{topo.colors}};
\code{\link{heat.colors}}; \code{\link{terrain.colors}}
}

\examples{
# Load the YFT diet data (of class diet)
data(yftdiet)  
class(yftdiet)

# Load the prey taxa data
data(PreyTaxonSort)
PreyTaxonSort

# Example where no palette is given
val <- apc(x = yftdiet, preyfile = PreyTaxonSort, palette = NULL, check = TRUE)
node.colsY <- val$cols
dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes
head(dietPP)

# Example where palette is given
pal <- c(topo.colors(10)[1:2], heat.colors(10)[1:2], terrain.colors(25)[1:8])
val <- apc(x = yftdiet, preyfile = PreyTaxonSort, palette = pal, check = TRUE)
node.colsY <- val$cols
dietPP <- val$x   # updated diet matrix with Group assigned prey taxa codes
head(dietPP)  


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~colors }
\keyword{ ~colours }% __ONLY ONE__ keyword per line
\keyword{ ~palette }% __ONLY ONE__ keyword per line
