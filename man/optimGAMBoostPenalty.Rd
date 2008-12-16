\name{optimGAMBoostPenalty}
\alias{optimGAMBoostPenalty}
\title{Coarse line search for adequate GAMBoost penalty parameter}
\description{
This routine helps in finding a penalty value that leads to an ``optimal'' number of boosting steps for GAMBoost (determined by AIC or cross-validation) that is not too small/in a specified range.  
}
\synopsis{
optimGAMBoostPenalty(x=NULL,y,x.linear=NULL,minstepno=50,maxstepno=200,start.penalty=500,method=c("AICmin","CVmin"),
                     penalty=100,penalty.linear=100,
                     just.penalty=FALSE,iter.max=10,upper.margin=0.05,trace=TRUE,parallel=FALSE,                               
                     calc.hat=TRUE,calc.se=TRUE,
                     which.penalty=ifelse(!is.null(x),"smoothness","linear"),\dots)
}
\usage{
optimGAMBoostPenalty(x=NULL,y,x.linear=NULL,
                     minstepno=50,maxstepno=200,method=c("AICmin","CVmin"),
                     which.penalty=ifelse(!is.null(x),"smoothness","linear"),
                     start.penalty=500,penalty=100,penalty.linear=100,
                     iter.max=10,upper.margin=0.05,
                     just.penalty=FALSE,trace=TRUE,parallel=FALSE,\dots)
}
\arguments{
\item{x}{\code{n * p} matrix of covariates with potentially non-linear influence. If this is not given (and argument \code{x.linear} is employed), a generalized linear model is fitted.}
\item{y}{response vector of length \code{n}.}
\item{x.linear}{optional \code{n * q} matrix of covariates with linear influence.}
\item{minstepno, maxstepno}{range of boosting steps in which the ``optimal'' number of boosting steps is wanted to be.}
\item{method}{determines how the optimal number of boosting steps corresponding to a fixed penalty is evaluated. With \code{"AICmin"} the AIC is used and with \code{"CVmin"} cross-validation is used as a criterion.}
\item{which.penalty}{indicates whether the penalty for the smooth components (value \code{"smoothness"}) or for the linear components (\code{"linear"}) should be optimized.}
\item{start.penalty}{start value for the search for the appropriate penalty.}
\item{penalty, penalty.linear}{penalty values for the smooth or parametric components, respectively for the type of components for which the penalty is not optimized.}
\item{just.penalty}{logical value indicating whether just the optimal penalty value should be returned or a \code{\link{GAMBoost}} fit performed with this penalty.}
\item{iter.max}{maximum number of search iterations.}
\item{upper.margin}{specifies the fraction of \code{maxstepno} which is used as an upper margin in which an AIC/cross-validation minimum is not taken to be one. (Necessary because of random fluctuations of these criteria).}
\item{trace}{logical value indicating whether information on progress should be printed.}
\item{parallel}{logical value indicating whether evaluation of cross-validation folds should be performed in parallel
on a compute cluster, when using \code{method="CVmin"}. This requires library \code{snowfall}.}
\item{\dots}{miscellaneous parameters for \code{\link{GAMBoost}}.}
}
\details{
The penalty parameter(s) for \code{\link{GAMBoost}} have to be chosen only very coarsely.  In Tutz and Binder (2006) it is suggested just to make sure, that the optimal number of boosting steps (according to AIC or cross-validation) is larger or equal to 50.  With a smaller number of steps boosting may become too ``greedy'' and show sub-optimal performance.  This procedure uses very a coarse line search and so one should specify a rather large range of boosting steps.  

Penalty optimization based on AIC should work fine most of the time, but for a large number of covariates (e.g. 500 with 100 observations) problems arise and (more costly) cross-validation should be employed.  
}
\value{
\code{GAMBoost} fit with the optimal penalty (with an additional component \code{optimGAMBoost.criterion}
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
\code{\link{GAMBoost}}
}
\examples{
\dontrun{
##  Generate some data 

x <- matrix(runif(100*8,min=-1,max=1),100,8)             
eta <- -0.5 + 2*x[,1] + 2*x[,3]^2
y <- rbinom(100,1,binomial()$linkinv(eta))

##  Find a penalty (starting from a large value, here: 5000) 
##  that leads to an optimal number of boosting steps (based in AIC) 
##  in the range [50,200] and return a GAMBoost fit with
##  this penalty

opt.gb1 <- optimGAMBoostPenalty(x,y,minstepno=50,maxstepno=200,
                                start.penalty=5000,family=binomial(),
                                trace=TRUE)

#   extract the penalty found/used for the fit
opt.gb1$penalty

}

}
\keyword{models} \keyword{smooth} \keyword{regression}
