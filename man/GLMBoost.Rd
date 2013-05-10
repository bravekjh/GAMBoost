\name{GLMBoost}
\alias{GLMBoost}
\title{Generalized linear model by likelihood based boosting}
\description{
\code{GLMBoost} a convenience wrapper around \code{\link{GAMBoost}}, for fitting generalized linear models by likelihood based boosting.
}
\usage{
GLMBoost(x,y,penalty=length(y),standardize=TRUE,\dots) 
}
\arguments{
\item{x}{\code{n * q} matrix of covariates with linear influence.}
\item{y}{response vector of length \code{n}.}
\item{penalty}{penalty value (scalar or vector of length q) for update of individual linear components in each boosting step. If this is set to \code{0} the covariates enter the model as mandatory covariates, which are updated together with the intercept term in each step.}
\item{standardize}{logical value indicating whether linear covariates should be standardized for estimation.}
\item{\dots}{arguments that should be passed to \code{\link{GAMBoost}}}
}
\value{
Object returned by call to \code{\link{GAMBoost}} (see documentation there), with additional class \code{GLMBoost}.
}
\author{
Harald Binder \email{binderh@uni-mainz.de}
}
\references{
Tutz, G. and Binder, H. (2007) Boosting ridge regression. \emph{Computational Statistics \& Data Analysis}, \bold{51}(12), 6044--6059.
}
\seealso{
\code{\link{GAMBoost}}, \code{\link{predict.GLMBoost}}.
}
\examples{
##  Generate some data 
x <- matrix(runif(100*8,min=-1,max=1),100,8)             
eta <- -0.5 + 2*x[,1] + 4*x[,3]
y <- rbinom(100,1,binomial()$linkinv(eta))

##  Fit a model with only linear components
gb1 <- GLMBoost(x,y,penalty=100,stepno=100,trace=TRUE,family=binomial()) 

#   Inspect the AIC for a minimum
plot(gb1$AIC) 

#   print the selected covariates, i.e., covariates with non-zero estimates
getGAMBoostSelected(gb1)

##  Make the first two covariates mandatory

gb2 <- GLMBoost(x,y,penalty=c(0,0,rep(100,ncol(x)-2)),
                stepno=100,family=binomial(),trace=TRUE) 

}
\keyword{models} \keyword{regression}
