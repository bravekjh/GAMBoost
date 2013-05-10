\name{cv.GLMBoost}
\alias{cv.GLMBoost}
\title{Cross-validation for GLMBoost fits}
\description{
Performs a convenience wrapper around \code{\link{cv.GAMBoost}} for performing a K-fold cross-validation for \code{\link{GLMBoost}} in search for the optimal number of boosting steps.  
}
\usage{
cv.GLMBoost(x,y,penalty=length(y),just.criterion=TRUE,...)
}
\arguments{
\item{y}{response vector of length \code{n}.}
\item{x}{\code{n * q} matrix of covariates with linear influence.}
\item{penalty}{penalty for the covariates with linear influence.}
\item{just.criterion}{logical value indicating wether a list with the goodness-of-fit information should be returned or a \code{GLMBoost} fit with the optimal number of steps.}
\item{\dots}{parameters to be passed to \code{\link{cv.GAMBoost}} or subsequently \code{\link{GAMBoost}}}
}
\value{
\code{GLMBoost} fit with the optimal number of boosting steps or list with the following components:
\item{criterion}{vector with goodness-of fit criterion for boosting step \code{1} , \dots , \code{maxstep}}
\item{se}{vector with standard error estimates for the goodness-of-fit criterion in each boosting step.}
\item{selected}{index of the optimal boosting step.}
}
\author{
Harald Binder \email{binderh@uni-mainz.de}
}
\seealso{
\code{\link{GLMBoost}}, \code{\link{cv.GAMBoost}}, \code{\link{GAMBoost}}
}
\examples{
\dontrun{
##  Generate some data 
x <- matrix(runif(100*8,min=-1,max=1),100,8)             
eta <- -0.5 + 2*x[,1] + 4*x[,3]
y <- rbinom(100,1,binomial()$linkinv(eta))

##  Fit the model with only linear components
gb1 <- GLMBoost(x,y,penalty=100,stepno=100,trace=TRUE,family=binomial()) 


##  10-fold cross-validation with prediction error as a criterion
gb1.crit <- cv.GLMBoost(x,y,penalty=100,maxstepno=100,trace=TRUE,
                        family=binomial(),
                        K=10,type="error")

##  Compare AIC and estimated prediction error

which.min(gb1$AIC)          
which.min(gb1.crit$criterion)
}

}
\keyword{models} \keyword{regression}
