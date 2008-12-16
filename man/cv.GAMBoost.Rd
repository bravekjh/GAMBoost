\name{cv.GAMBoost}
\alias{cv.GAMBoost}
\title{Cross-validation for GAMBoost fits}
\description{
Performs a K-fold cross-validation for \code{\link{GAMBoost}} in search for the optimal number of boosting steps.  
}
\synopsis{
cv.GAMBoost(x=NULL,y,x.linear=NULL,maxstepno=500,family=binomial(),weights=rep(1,length(y)),
            calc.hat=TRUE,calc.se=TRUE,trace=FALSE,parallel=FALSE,upload.x=TRUE,
            K=10,type=c("loglik","error","L2"),pred.cutoff=0.5,just.criterion=FALSE,\dots) 
}
\usage{
cv.GAMBoost(x=NULL,y,x.linear=NULL,maxstepno=500,
            K=10,type=c("loglik","error","L2"),pred.cutoff=0.5,
            just.criterion=FALSE,trace=FALSE,parallel=FALSE,upload.x=TRUE,\dots) 
}
\arguments{
\item{x}{\code{n * p} matrix of covariates with potentially non-linear influence. If this is not given (and argument \code{x.linear} is employed), a generalized linear model is fitted.}
\item{y}{response vector of length \code{n}.}
\item{x.linear}{optional \code{n * q} matrix of covariates with linear influence.}
\item{maxstepno}{maximum number of boosting steps to evaluate.}
\item{K}{number of folds to be used for cross-validation.}
\item{type, pred.cutoff}{goodness-of-fit criterion: likelihood (\code{"loglik"}), error rate for binary response data (\code{"error"}) or squared error for others (\code{"L2"}).
For binary response data and the \code{"error"} criterion \code{pred.cutoff} specifies the p value cutoff for prediction of class 1 vs 0.}
\item{just.criterion}{logical value indicating wether a list with the goodness-of-fit information should be returned or a \code{GAMBoost} fit with the optimal number of steps.}
\item{trace}{logical value indicating whether information on progress should be printed.}
\item{parallel}{logical value indicating whether evaluation of cross-validation folds should be performed in parallel
on a compute cluster. This requires library \code{snowfall}.}
\item{upload.x}{logical value indicating whether \code{x} and \code{x.linear} should/have to be uploaded to the
compute cluster for parallel computation. Uploading these only once (using \code{sfExport(x,x.linear)} from library \code{snowfall}) can save much time for large data sets.}
\item{\dots}{miscellaneous parameters for the calls to \code{\link{GAMBoost}}}
}
\value{
\code{GAMBoost} fit with the optimal number of boosting steps or list with the following components:
\item{criterion}{vector with goodness-of fit criterion for boosting step \code{1} , \dots , \code{maxstep}}
\item{se}{vector with standard error estimates for the goodness-of-fit criterion in each boosting step.}
\item{selected}{index of the optimal boosting step.}
}
\author{
Harald Binder \email{binderh@fdm.uni-freiburg.de}
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

##  Fit the model with smooth components

gb1 <- GAMBoost(x,y,penalty=400,stepno=100,trace=TRUE,family=binomial()) 

##  10-fold cross-validation with prediction error as a criterion

gb1.crit <- cv.GAMBoost(x,y,penalty=400,maxstepno=100,trace=TRUE,
                        family=binomial(),
                        K=10,type="error",just.criterion=TRUE)

##  Compare AIC and estimated prediction error

which.min(gb1$AIC)          
which.min(gb1.crit$criterion)
}

}
\keyword{models} \keyword{smooth} \keyword{regression}
