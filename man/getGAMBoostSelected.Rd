\name{getGAMBoostSelected}
\alias{getGAMBoostSelected}
\title{Identify selected/significant covariates from a GAMBoost object}
\description{
Extracts the information from a GAMBoost object which covariates have received any update up to a specific boosting step and for which smooth estimates the pointwise confidence bands do not contain the zero line.  
}
\usage{
getGAMBoostSelected(object,at.step=NULL)
}
\arguments{
\item{object}{fitted GAMBoost object from a \code{\link{GAMBoost}} call.}
\item{at.step}{boosting step for which the information should be extracted. If none is given, the final boosting step is examined.}
}
\value{
\item{smooth}{indices of smooth components which received any update up to the given step.}
\item{smoothbands}{indices of smooth components for which the pointwise confidence bands at the given step do not contain the null line (only available when confidence bands have been fitted; see option \code{se.fit} of \code{\link{GAMBoost}}).}
\item{parametric}{indices of parametric components which received any update up to the given step.}
}
\author{
Harald Binder \email{binderh@uni-mainz.de}
}
\seealso{
\code{\link{GAMBoost}}, \code{\link{plot.GAMBoost}}
}
\examples{
## see examples for 'GAMBoost' and 'plot.GAMBoost'
}
\keyword{models} \keyword{smooth} \keyword{regression}
