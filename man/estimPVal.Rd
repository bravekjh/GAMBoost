\name{estimPVal}
\alias{estimPVal}
\title{Estimate p-values for a model fitted by GAMBoost or GLMBoost}
\description{
Performs permutation-based p-value estimation for the optional covariates in a fit from \code{\link{GAMBoost}} or \code{\link{GAMBoost}}. Currently binary response models with linear effects are supported, and the components have to be selected with \code{criterion="score"}}
\usage{
estimPVal(object,x,y,permute.n=10,per.covariate=FALSE,parallel=FALSE,
          multicore=FALSE,trace=FALSE,...)
}
\arguments{
\item{object}{fit object obtained from \code{\link{GAMBoost}} or \code{\link{GLMBoost}}.}
\item{x}{\code{n * p} matrix of covariates with linear effect. This has to be the same that was used as \code{x.linear} in the call to \code{\link{GAMBoost}} or \code{x} in \code{\link{GLMBoost}}.}
\item{y}{response vector. This has to be the same that was used in the call to \code{\link{GAMBoost}} or \code{\link{GLMBoost}}.}
\item{permute.n}{number of permutations employed for obtaining a null distribution.}
\item{per.covariate}{logical value indicating whether a separate null distribution should be considered for each covariate. A larger number of permutations will be needed if this is wanted.}
\item{parallel}{logical value indicating whether computations for obtaining a null distribution via permutation should be performed in parallel on a compute cluster. Parallelization is performed via the package \code{snowfall} and the initialization function of of this package, \code{sfInit}, should be called before calling \code{estimPVal}.}
\item{multicore}{indicates whether computations in the permuted data sets should be performed in parallel, using package \code{multicore}. If \code{TRUE}, package \code{multicore} is employed using the default number of cores. A value larger than \code{1} is taken to be the number of cores that should be employed.}
\item{trace}{logical value indicating whether progress in estimation should be indicated by printing the number of the permutation that is currently being evaluated.}
\item{\dots}{miscellaneous parameters for the calls to \code{\link{GAMBoost}}}
}
\details{
As p-value estimates are based on permutations, random numbers are drawn for determining permutation indices. Therfore, the results depend on the state of the random number generator. This can be used to explore the variability due to random variation and help to determine an adequate value for \code{permute.n}. A value of 100 should be sufficient, but this can be quite slow. If there is a considerable number of covariates, e.g., larger than 100, a much smaller number of permutations, e.g., 10, might already work well. The estimates might also be negatively affected, if only a small number of boosting steps (say <50) was employed for the original fit. 
}
\value{
Vector with p-value estimates, one value for each optional covariate with linear effect specificed in the original call to \code{\link{GAMBoost}} or \code{\link{GLMBoost}}.
}
\author{
Harald Binder \email{binderh@uni-mainz.de}
}
\references{
Binder, H., Porzelius, C. and Schumacher, M. (2009). Rank-based p-values for sparse high-dimensional risk prediction models fitted by componentwise boosting. FDM-Preprint Nr. 101, University of Freiburg, Germany.
}
\seealso{
\code{\link{GAMBoost}}, \code{\link{GLMBoost}}
}
\examples{
\dontrun{
##  Generate some data 
x <- matrix(runif(100*8,min=-1,max=1),100,8)             
eta <- -0.5 + 2*x[,1] + 4*x[,3]
y <- rbinom(100,1,binomial()$linkinv(eta))

##  Fit a model with only linear components
gb1 <- GLMBoost(x,y,penalty=100,stepno=100,trace=TRUE,family=binomial(),criterion="score") 

#   estimate p-values

p1 <- estimPVal(gb1,x,y,permute.n=10)

#   get a second vector of estimates for checking how large
#   random variation is

p2 <- estimPVal(gb1,x,y,permute.n=10)

plot(p1,p2,xlim=c(0,1),ylim=c(0,1),xlab="permute 1",ylab="permute 2")
}

}
\keyword{models} \keyword{regression} \keyword{survial}
