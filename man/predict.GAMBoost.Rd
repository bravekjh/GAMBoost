\name{predict.GAMBoost}
\alias{predict.GAMBoost}
\title{Predict method for GAMBoost fits}
\description{
Obtains predictions at specified boosting steps from a GAMBoost object fitted by \code{\link{GAMBoost}}.  
}
\usage{
\method{predict}{GAMBoost}(object,newdata=NULL,newdata.linear=NULL,
                 at.step=NULL,type=c("link","response","terms"),\dots)
}
\arguments{
\item{object}{fitted GAMBoost object from a \code{\link{GAMBoost}} call.}
\item{newdata}{\code{n.new * p} matrix with new covariate values for smooth components. If just prediction for the training data is wanted or just a generalized linear model has been fitted, it can be omitted.}
\item{newdata.linear}{matrix with new covariate values for linear components. If linear components have been fitted and this is not given, the contribution of the linear components will be ignored for prediction.}
\item{at.step}{scalar or vector of boosting step(s) at which prediction is wanted. If \code{type="terms"} is used, only one step is admissible. If no step is given, the final boosting step is used.}
\item{type}{type of prediction to be returned: \code{"link"} gives prediction at the level of the predictor, \code{"response"} at the response level. \code{"terms"} returns individual contributions of the smooth components to the predictor.}
\item{\dots}{miscellaneous arguments, none of which is used at the moment.}
}
\value{
For \code{type="link"} and \code{type="response"} a vector of length \code{n.new} (\code{at.step} being a scalar) or a \code{n.new * length(at.step)} matrix (\code{at.step} being a vector) with predictions is returned.  
For \code{type="terms"} a \code{n.new * p+1} matrix with contributions of the smooth components to the predictor is returned.  
}
\author{
Harald Binder \email{binderh@fdm.uni-freiburg.de}
}
\examples{
##  Generate some data 
x <- matrix(runif(100*3,min=-1,max=1),100,3)             
eta <- -0.5 + 2*x[,1] + 4*x[,3]^2
y <- rbinom(100,1,binomial()$linkinv(eta))

##  Fit the model with smooth components
gb1 <- GAMBoost(x,y,penalty=200,stepno=100,trace=TRUE,family=binomial()) 

##  Extract predictions

#   at final boosting step
predict(gb1,type="response")

#   at 'optimal' boosting step (with respect to AIC)
predict(gb1,at.step=which.min(gb1$AIC),type="response")

#   matrix with predictions at predictor level for all boosting steps
predict(gb1,at.step=1:100,type="link")

}
\keyword{models} \keyword{smooth} \keyword{regression}
