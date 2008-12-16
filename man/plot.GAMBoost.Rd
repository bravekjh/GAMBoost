\name{plot.GAMBoost}
\alias{plot.GAMBoost}
\title{Plots of the smooth functions from a GAMBoost fit}
\description{
Generates plots for the smooth components from a \code{\link{GAMBoost}} fit at a specific boosting step.
}
\synopsis{
\method{plot}{GAMBoost}(x,select=NULL,at.step=NULL,add=FALSE,phi=1,ylim=NULL,xlab=NULL,ylab=NULL,\dots)
}
\usage{
\method{plot}{GAMBoost}(x,select=NULL,at.step=NULL,add=FALSE,phi=1,\dots)
}
\arguments{
\item{x}{fitted GAMBoost object from a \code{\link{GAMBoost}} call.}
\item{select}{indices of the smooth component(s) for which plots should be generated. If none are specified, all are used.}
\item{at.step}{boosting step from which the estimates for the smooth functions should be evaluated. If not given, the final boosting step is used.}
\item{add}{logical value indicating whether the plot(s) should be added to the current plot.}
\item{phi}{scale parameter for the confidence bands.}
\item{\dots}{miscellaneous plotting parameters given to the low level plotting routine.}
}
\value{
A plot is produced for the specified smooth components in the \code{GAMBoost} fit.  Pointwise confidence bands are plotted when the standard error information has been calculated (option \code{calc.se=TRUE} in the call to \code{GAMBoost}).  
}
\author{
Harald Binder \email{binderh@fdm.uni-freiburg.de}
}
\examples{
##  Generate some data 

x <- matrix(runif(100*8,min=-1,max=1),100,8)             
eta <- -0.5 + 2*x[,1] + 2*x[,3]^2
y <- rbinom(100,1,binomial()$linkinv(eta))

##  Fit the model with smooth components
gb1 <- GAMBoost(x,y,penalty=400,stepno=100,trace=TRUE,family=binomial()) 

##  Plot smooth components of fit

#   all, at final boosting step
par(mfrow=c(2,4))
plot(gb1)

#   components that received an update up to the 'optimal' boosting step

selected <- getGAMBoostSelected(gb1,at.step=which.min(gb1$AIC))

par(mfrow=c(1,length(selected$smooth)))
plot(gb1,select=selected$smooth)

#   components where the estimate at the 'optimal' boosting step does not
#   contain the null line
par(mfrow=c(1,length(selected$smoothbands)))
plot(gb1,select=selected$smoothbands)

}
\keyword{models} \keyword{smooth} \keyword{regression}
