\name{igam}
\alias{igam}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Isotope Generalized Additive Model a_5 b^1
}

\description{
Fits a Generalized Additive Model (GAM). This function is a wrapper for the
\code{mgcv} package developed by Simon Wood.
}
\usage{
igam(formula, data, control = list(), G = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A GAM formula (see \code{formula.gam} and also \code{gam.models}). This is exactly 
like the formula for a GLM except that smooth terms, \code{s} and \code{te} 
can be added to the right hand side to specify that the linear predictor 
depends on smooth functions of predictors (or linear functionals of these).
}
  \item{data}{
A data frame or list containing the model response variable and covariates 
required by the formula. By default the variables are taken from environment
(formula): typically the environment from which gam is called.
}
  \item{control}{
A list of fit control parameters to replace defaults returned by gam.control. 
Values not set assume default values.
}
  \item{G}{
Usually NULL, but may contain the object returned by a previous call to gam with 
fit=FALSE, in which case all other arguments are ignored except for gamma, 
in.out, scale, control, method optimizer and fit.
}
  \item{\dots}{
Describe \code{\dots} here
}
}
\details{
This function is a wrapper for Simon Wood's \code{gam} function in the \code{mgcv}
package. For specific details regarding GAMs, see \code{\link{mgcv}}. By default,
the family is assumed Gaussian. If an alternative family is required, see the
\code{mgcv} package and refit the model using the \code{gam} function.
}
\value{
See \code{igam.object} for more details.
}
\references{
Kuhnert, P.M., Duffy, L. M and Olson, R.J. (2012) The Analysis of Predator Diet 
and Stable Isotope Data, Journal of Statistical Software, In Prep. 

Wood, S.N. (2011) Fast stable restricted maximum likelihood and marginal 
likelihood estimation of semiparametric generalized linear models. 
Journal of the Royal Statistical Society (B) 73(1):3-36

Wood, S.N. (2004) Stable and efficient multiple smoothing parameter estimation 
for generalized additive models. J. Amer. Statist. Ass. 99:673-686.

Wood, S.N. (2006) Generalized Additive Models: an introduction with R, CRC
}




\seealso{
\code{\link{mgcv}}
}
\examples{
# Load Isotope Dataset
#data(yftISO)

# Fit GAM
#Yd15N.gam <- igam(d15N ~ Year + s(FL_mm) + s(SST) + s(Lon, Lat),
#                     data = yftISO)
#summary(Yd15N.gam)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~gam }
\keyword{ ~isotope }% __ONLY ONE__ keyword per line
