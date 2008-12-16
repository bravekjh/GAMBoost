\name{optimGLMBoostPenalty}
\alias{optimGLMBoostPenalty}
\title{Coarse line search for adequate GLMBoost penalty parameter}
\description{
This routine is a convenience wrapper around \code{\link{optimGAMBoostPenalty}} for finding a penalty value that leads to an ``optimal'' number of boosting steps for GLMBoost (determined by AIC or cross-validation) that is not too small/in a specified range.  
}
\usage{
optimGLMBoostPenalty(x,y,start.penalty=length(y),just.penalty=FALSE,...)}
\arguments{
\item{x}{\code{n * q} matrix of covariates with linear influence.}
\item{y}{response vector of length \code{n}.}
\item{start.penalty}{start value for the search for the appropriate penalty.}
\item{just.penalty}{logical value indicating whether just the optimal penalty value should be returned or a \code{\link{GLMBoost}} fit performed with this penalty.}
\item{\dots}{miscellaneous parameters for \code{\link{optimGAMBoostPenalty}}.}
}
\details{
The penalty parameter(s) for \code{\link{GLMBoost}} have to be chosen only very coarsely.  In Tutz and Binder (2006) it is suggested just to make sure, that the optimal number of boosting steps (according to AIC or cross-validation) is larger or equal to 50.  With a smaller number of steps boosting may become too ``greedy'' and show sub-optimal performance.  This procedure uses very a coarse line search and so one should specify a rather large range of boosting steps.  

Penalty optimization based on AIC should work fine most of the time, but for a large number of covariates (e.g. 500 with 100 observations) problems arise and (more costly) cross-validation should be employed.  
}
\value{
\code{GLMBoost} fit with the optimal penalty (with an additional component \code{optimGAMBoost.criterion}
giving the values of the criterion (AIC or cross-validation) corresponding to the final penalty) 
or just the optimal penalty value itself.  
}
\author{
Written by Harald Binder \email{binderh@fdm.uni-freiburg.de}, matching closely the original Fortran implementation employed
for Tutz and Binder (2006). 
}
\references{
Tutz, G. and Binder, H. (2006) Generalized additive modelling with implicit variable selection by likelihood based boosting. \emph{Biometrics}, \bold{51}, 961--971.
}
\seealso{
\code{\link{GLMBoost}}, \code{\link{optimGAMBoostPenalty}}, \code{\link{GAMBoost}}
}
\examples{
\dontrun{
##  Generate some data 

##  Generate some data 
x <- matrix(runif(100*8,min=-1,max=1),100,8)             
eta <- -0.5 + 2*x[,1] + 4*x[,3]
y <- rbinom(100,1,binomial()$linkinv(eta))

##  Find a penalty (starting from a large value, here: 5000) 
##  that leads to an optimal number of boosting steps (based in AIC) 
##  in the range [50,200] and return a GLMBoost fit with
##  this penalty

opt.gb1 <- optimGLMBoostPenalty(x,y,minstepno=50,maxstepno=200,
                                start.penalty=5000,family=binomial(),
                                trace=TRUE)

#   extract the penalty found/used for the fit
opt.gb1$penalty

}

}
\keyword{models} \keyword{regression}
