\name{optimStepSizeFactor}
\alias{optimStepSizeFactor}
\title{Coarse line search for optimum step-size modification factor}
\description{
This routine helps in finding an optimum step-size modification factor for \code{\link{GAMBoost}}, i.e., that results in an optimum in terms of cross-validated log-likelihood.  
}
\usage{
optimStepSizeFactor(x=NULL,y,x.linear=NULL,
                    direction=c("down","up","both"),start.stepsize=0.1,
                    iter.max=10,constant.cv.res=NULL,parallel=FALSE,
                    trace=FALSE,...) 
}
\arguments{
\item{x}{\code{n * p} matrix of covariates with potentially non-linear influence. If this is not given (and argument \code{x.linear} is employed), a generalized linear model is fitted.}
\item{y}{response vector of length \code{n}.}
\item{x.linear}{optional \code{n * q} matrix of covariates with linear influence.}
\item{direction}{direction of line search for an optimal step-size modification factor (starting from value 1).}
\item{start.stepsize}{step size used for the line search. A final step is performed using half this size.}
\item{iter.max}{maximum number of search iterations.}
\item{constant.cv.res}{result of \code{\link{cv.GAMBoost}} (with \code{just.criterion=TRUE}) for \code{stepsize.factor.linear=1}, that can be provided for saving computing time, if it already is available.}
\item{parallel}{logical value indicating whether evaluation of cross-validation folds should be performed in parallel
on a compute cluster. This requires library \code{snowfall}.}
\item{trace}{logical value indicating whether information on progress should be printed.}
\item{\dots}{miscellaneous parameters for \code{\link{cv.GAMBoost}}.}
}
\details{
A coarse line search is performed for finding the best parameter \code{stepsize.factor.linear} for \code{\link{GAMBoost}}. If an \code{pendistmat.linear} argument is provided (which is passed on to \code{\link{GAMBoost}}), a search for factors smaller than 1 is sensible (corresponding to \code{direction="down"}). If no connection information is provided, it is reasonable to employ \code{direction="both"}, for avoiding restrictions without subject matter knowledge.
}
\value{
List with the following components:
\item{factor.list}{array with the evaluated step-size modification factors.}
\item{critmat}{matrix with the mean log-likelihood for each step-size modification factor in the course of the boosting steps.}
\item{optimal.factor.index}{index of the optimal step-size modification factor.}
\item{optimal.factor}{optimal step-size modification factor.}
\item{optimal.step}{optimal boosting step number, i.e., with minimum mean log-likelihood, for step-size modification factor \code{optimal.factor}.}
}
\author{
Written by Harald Binder \email{binderh@fdm.uni-freiburg.de}. 
}
\references{
Binder, H. and Schumacher, M. (2009). Incorporating pathway information into boosting estimation of high-dimensional risk prediction models. BMC Bioinformatics. 10:18.
}
\seealso{
\code{\link{GAMBoost}}, \code{\link{cv.GAMBoost}}
}
\examples{
\dontrun{
##  Generate some data 
n <- 100; p <- 10

#   covariates with non-linear (smooth) effects
x <- matrix(runif(n*p,min=-1,max=1),n,p)             
eta <- -0.5 + 2*x[,1] + 2*x[,3]^2 + x[,9]-.5
y <- rbinom(n,1,binomial()$linkinv(eta))

#  Determine step-size modification factor for a generalize linear model
#  As there is no connection matrix,  perform search into both directions 

optim.res <- optimStepSizeFactor(direction="both",
                                y=y,x.linear=x,family=binomial(),
                                penalty.linear=200,
                                trace=TRUE)

#   Fit with obtained step-size modification parameter and optimal number of boosting
#   steps obtained by cross-validation

gb1 <- GAMBoost(x=NULL,y=y,x.linear=x,family=binomial(),penalty.linear=200,
                stepno=optim.res$optimal.step,
                stepsize.factor.linear=optim.res$optimal.factor) 

summary(gb1)

}

}
\keyword{models} \keyword{smooth} \keyword{regression}
