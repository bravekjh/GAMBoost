\name{predict.GLMBoost}
\alias{predict.GLMBoost}
\title{Predict method for GLMBoost fits}
\description{
Convienience wrapper for \code{\link{predict.GAMBoost}}, for obtaining predictions at specified boosting steps from a GAMBoost object fitted by \code{\link{GLMBoost}}.  
}
\usage{
\method{predict}{GLMBoost}(object,newdata=NULL,\dots)
}
\arguments{
\item{object}{fitted GAMBoost object from a \code{\link{GLMBoost}} call.}
\item{newdata}{\code{n.new * q} matrix with new covariate values for linear components. If just prediction for the training data is wanted, it can be omitted.}
\item{\dots}{arguments that should be passed to \code{\link{predict.GAMBoost}}.}
}
\value{
Value returned by \code{\link{predict.GAMBoost}} (see documentation there).
}
\author{
Harald Binder \email{binderh@fdm.uni-freiburg.de}
}
\examples{
##  Generate some data 
x <- matrix(runif(100*8,min=-1,max=1),100,8)             
eta <- -0.5 + 2*x[,1] + 4*x[,3]
y <- rbinom(100,1,binomial()$linkinv(eta))

##  Fit the model with only linear components
gb1 <- GLMBoost(x,y,penalty=100,stepno=100,trace=TRUE,family=binomial()) 

##  Extract predictions

#   at final boosting step
predict(gb1,type="response")

#   at 'optimal' boosting step (with respect to AIC)
predict(gb1,at.step=which.min(gb1$AIC),type="response")

#   matrix with predictions at predictor level for all boosting steps
predict(gb1,at.step=1:100,type="link")
}
\seealso{
\code{\link{GLMBoost}}, \code{\link{GAMBoost}}, \code{\link{predict.GAMBoost}}.
}
\keyword{models} \keyword{regression}
