\name{GAMBoost}
\alias{GAMBoost}
\title{Generalized additive model by likelihood based boosting}
\description{
\code{GAMBoost} is used to fit a generalized additive model by likelihood based boosting.  
It is especially suited for models with a large number of predictors with potentially 
non-linear influence. It provides smooth function estimates of covariate influence functions
together with confidence bands and approximate degrees of freedom.  
}
\usage{
GAMBoost(x=NULL,y,xmin=NULL,xmax=NULL,penalty=100,bdeg=2,pdiff=1,
         x.linear=NULL,standardize.linear=TRUE,penalty.linear=0,
         weights=rep(1,length(y)),stepno=500,family=binomial(),
         sparse.boost=FALSE,sparse.weight=1,calc.hat=TRUE,calc.se=TRUE,
         AIC.type="corrected",trace=FALSE) 
}
\arguments{
\item{x}{\code{n * p} matrix of covariates with potentially non-linear influence. If this is not given (and argument \code{x.linear} is employed), a generalized linear model is fitted.}
\item{y}{response vector of length \code{n}.}
\item{xmin, xmax}{optional vectors of length \code{p} specifying the lower and upper bound for the range of the smooth functions to be fitted.}
\item{penalty}{penalty value for the update of an individual smooth function in each boosting step.}
\item{bdeg, pdiff}{degree of the B-spline basis to be used for fitting smooth functions and difference of the
coefficient estimates to which the penalty should be applied.}
\item{x.linear}{optional \code{n * q} matrix of covariates with linear influence.}
\item{standardize.linear}{logical value indicating whether linear covariates should be standardized for estimation.}
\item{penalty.linear}{penalty value (scalar or vector of length q) for update of individual linear components in each boosting step. If this is set to \code{0} the covariates in \code{x.linear} enter the model as mandatory covariates, which are updated together with the intercept term in each step.}
\item{weights}{an optional vector of weights to be used in the fitting process.}
\item{stepno}{number of boosting steps (\code{m}).}
\item{family}{a description of the error distribution to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) Note that \code{GAMBoost} supports only canonical link functions and no scale parameter estimation, so \code{gaussian()}, \code{binomial()} and \code{poisson()} are the most plausible candidates.}
\item{sparse.boost}{logical value indicating whether a criterion considering degrees of freedom (specifically AIC) should be used for selecting a covariate for an update in each boosting step (sparse boosting), instead of the deviance.}
\item{sparse.weight}{factor modifying how the degrees of freedom enter into calculation of AIC for sparse boosting.}
\item{calc.hat}{logical value indicating whether the hat matrix should be computed for each boosting step. If set to \code{FALSE} no degrees of freedom and therefore e.g. no AIC will be available. On the other hand fitting will be faster (especially for a large number of observations).}
\item{calc.se}{logical value indicating whether confidence bands should be calculated. Switch this of for faster fitting.}
\item{AIC.type}{type of model selection criterion to be calculated (and also to be used for sparse boosting if applicable): in the Gaussian case \code{"classical"} gives AIC and \code{"corrected"} results in corrected AIC (Hurvich, Simonoff and Tsai, 1998); for all other response families standard AIC is used.}
\item{trace}{logical value indicating whether progress in estimation should be indicated by printing the index  of the covariate updated (\code{1}, \dots , \code{p} for smooth components and \code{p+1}, \dots ,\code{p-q}} for parametric components) in the current boosting step.}
}
\details{
The idea of likelihood based boosting (Tutz and Binder, 2006) is most easily understood for models with a Gaussian response.  There it results in repeated fitting of residuals (This idea is then transferred to the generalized case).  For obtaining an additive model GAMBoost uses a large number of boosting steps where in each step a penalized B-spline (of degree \code{bdeg}) (Eilers and Marx, 1996) is fitted to one covariate, the response being the residuals from the last step.  The covariate to be updated is selected by deviance (or in case of sparse boosting by some model selection criterion).  The B-spline coefficient estimates in each step are fitted under the constraint of a large penalty on their first (or higher order) differences. So in each step only a small adjustment is made.  Summing over all steps for each covariate a smooth function estimate is obtained. When no basis expansion is used, i.e. just coefficients of covariates are updated, (generalized) linear models are obtained.  

The main parameter of the algorithm is the number of boosting steps, given that the penalty is chosen large enough (If too small the minimum AIC will occur for a very early boosting step; see \code{\link{optimGAMBoostPenalty}}).  When there is a large number of covariates with potentially non-linear effect, having a single parameter (with adaptive smoothness assignment to single components performed automatically by the algorithm) is a huge advantage compared to approaches, where a smoothing parameter has to be selected for each single component. The biggest advantage over conventional methods for fitting generalized additive models (e.g. \code{mgcv:gam} or \code{gam:gam}) will therefore be obtained for a large number of covariates compared to the number of observations (e.g. 10 covariates with only 100 observations).  In addition GAMBoost performs well compared to other approaches when there is a small signal-to-noise ratio and/or a response (e.g. binary) with a small amount of information.

Note that the degrees of freedom (and based on these AIC and BIC) are just approximations, which are completely valid for example only when the order and the indices of the components updated is fixed.  This leads to problems especially when there is a very large number of covariates (e.g. 10000 covariates with only 100 observations).  Then it might be better (but also slower) to rely on cross validation (see \code{\link{cv.GAMBoost}}) for selection of the number of boosting steps.  

Note that the \code{gamboost} routine in the R package \code{mboost} implements a different kind of boosting strategy: gradient based boosting instead of likelihood based boosting.  The two approaches coincide only in special cases (e.g. L2 loss  and Gaussian response).  While gradient based boosting is more general, only likelihood based boosting allows e.g. for easily obtainable pointwise confidence bands.  
}
\value{
\code{GAMBoost} returns an object of class \code{GAMBoost}.  \code{GAMBoost} objects can be examined by print, summary, and plot.  \code{\link{getGAMBoostSelected}} can be used on them to identify selected/significant covariates.  

\item{x, x.linear, n}{original covariate values and number of observations.}
\item{penalty, penalty.linear}{penalties used in updating smooth and linear components.}
\item{stepno}{number of boosting steps.}
\item{family}{response family.}
\item{AIC.type}{type of AIC given in component \code{AIC} (applies only for Gaussian response).}
\item{deviance, trace, AIC, BIC}{vectors of length \code{m} giving deviance, approximate degrees of freedom, AIC, and BIC for each boosting step.}
\item{selected}{vector of length \code{m} given the index of the covariate updated in each boosting step (\code{1-p} for smooth components and \code{(p+1)}, ... , \code{(p-q)} for parametric components).}
\item{beta}{list of length \code{p+1}, each element being a matrix with \code{m+1} rows giving the estimated coefficients for the intercept term (\code{beta[[1]]}) and the smooth terms (\code{beta[[2]]}, ... , \code{beta[[p+1]]}).}
\item{beta.linear}{\code{m * q} matrix containing the coefficient estimates for the (standardized) linear covariates.}
\item{mean.liner, sd.linear}{vector of mean values and standard deviations used for standardizing the linear covariates.}
\item{hatmatrix}{hat matrix at the final boosting step.}
\item{eta}{\code{n * (m+1)} matrix with predicted value (at predictor level) for each boosting step.}
\item{predictors}{list of length \code{p+1} containing information (as a list) on basis expansions for the smooth model components (\code{predictors[[2]]}, ... , \code{predictors[[p+1]]}).}
}
\author{
Written by Harald Binder \email{binderh@fdm.uni-freiburg.de}, matching closely the original Fortran implementation employed
for Tutz and Binder (2006). 
}
\references{
Hurvich, C. M., Simonoff, J. S. and Tsai, C. L. (1998). Smoothing parameter selection in nonparametric regression using and improved Akaike information criterion. Journal of the Royal Statistical Society B, \bold{60}(2), 271--293.

Eilers, P. H. C. and Marx, B. D. (1996) Flexible smoothing with B-splines and penalties. \emph{Statistical Science}, \bold{11}(2), 89--121.

Tutz, G. and Binder, H. (2006) Generalized additive modelling with implicit variable selection by likelihood based boosting. \emph{Biometrics}. To appear.

Tutz, G. and Binder, H. (2005) Boosting ridge regression. Discussion Paper 418, SFB 386, Ludwig-Maximilians-University Munich.
}
\seealso{
\code{\link{getGAMBoostSelected}}, \code{\link{plot.GAMBoost}}, \code{\link{predict.GAMBoost}}, \code{\link{optimGAMBoostPenalty}}.
}
\examples{
##  Generate some data 
n <- 100; p <- 8; q <- 2

#   covariates with non-linear (smooth) effects
x <- matrix(runif(n*p,min=-1,max=1),n,p)             

#   binary covariates
x.linear <- matrix(round(runif(n*q,min=0,max=1)),n,q)

#   1st and 3rd smooth covariate and 1st linear covariate are informative
eta <- -0.5 + 2*x[,1] + 2*x[,3]^2 + x.linear[,1]-.5

y <- rbinom(n,1,binomial()$linkinv(eta))

##  Fit a model with just smooth components
gb1 <- GAMBoost(x,y,penalty=500,stepno=100,family=binomial(),trace=TRUE) 

#   Inspect the AIC for a minimum
plot(gb1$AIC) # still falling at boosting step 100 so we need more steps
              # or a smaller penalty (use 'optimGAMBoostPenalty' for
              # automatic penalty optimization)

##  Include two binary covariates as mandatory without penalty
##  (appropriate for example for 'treatment/control')
##  modelled as 'linear' predictors

gb2 <- GAMBoost(x,y,penalty=200,
                x.linear=x.linear,penalty.linear=0,
                stepno=100,family=binomial(),trace=TRUE) 

##  Include first binary covariates as mandatory and second
##  as optional (e.g 'treatment/control' and 'female/male')

gb3 <- GAMBoost(x,y,penalty=200,
                x.linear=x.linear,penalty.linear=c(0,100),
                stepno=100,family=binomial(),trace=TRUE) 

#   Get summary with fitted covariates and estimates for
#   the parametric components
summary(gb3)

#   Extract boosted components at 'optimal' boosting step
selected <- getGAMBoostSelected(gb3,at.step=which.min(gb3$AIC))

#   Plot all smooth components at final boosting step
par(mfrow=c(2,4))
plot(gb3)

#   plot smooth components for which the null line is not inside the bands
#   at 'optimal' boosting step (determined by AIC)
par(mfrow=c(1,length(selected$smoothbands)))
plot(gb3,select=selected$smoothbands,at.step=which.min(gb3$AIC))

##   Fit a generalized linear model for comparison

x.linear <- cbind(x,x.linear)
gb4 <- GAMBoost(x=NULL,y=y,x.linear=x.linear,penalty.linear=100,
                stepno=100,trace=TRUE,family=binomial())

#   Compare with generalized additive model fit
plot(gb3$AIC,type="l"); lines(1:100,gb4$AIC,lty=2)

}
\keyword{models} \keyword{smooth} \keyword{regression}
