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
         x.linear=NULL,standardize.linear=TRUE,
         penalty.linear=0,subset=NULL,
         criterion=c("deviance","score"),stepsize.factor.linear=1,
         sf.scheme=c("sigmoid","linear"),pendistmat.linear=NULL,
         connected.index.linear=NULL,
         weights=rep(1,length(y)),stepno=500,family=binomial(),
         sparse.boost=FALSE,sparse.weight=1,calc.hat=TRUE,calc.se=TRUE,
         AIC.type=c("corrected","classical"),return.score=TRUE,trace=FALSE) 
}
\arguments{
\item{x}{\code{n * p} matrix of covariates with potentially non-linear influence. If this is not given (and argument \code{x.linear} is employed), a generalized linear model is fitted.}
\item{y}{response vector of length \code{n}.}
\item{xmin, xmax}{optional vectors of length \code{p} specifying the lower and upper bound for the range of the smooth functions to be fitted.}
\item{penalty}{penalty value for the update of an individual smooth function in each boosting step.}
\item{bdeg, pdiff}{degree of the B-spline basis to be used for fitting smooth functions and difference of the
coefficient estimates to which the penalty should be applied. When \code{pdiff} is a vector, the penalties corresponding to single vector elements are combined to enforce several types of smoothness simultaneously.}
\item{x.linear}{optional \code{n * q} matrix of covariates with linear influence.}
\item{standardize.linear}{logical value indicating whether linear covariates should be standardized for estimation.}
\item{penalty.linear}{penalty value (scalar or vector of length q) for update of individual linear components in each boosting step. If this is set to \code{0} the covariates in \code{x.linear} enter the model as mandatory covariates, which are updated together with the intercept term in each step.}
\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
\item{criterion}{Indicates the criterion to be employed for selecting the best update in each boosting step for the linear components. \code{deviance} corresponds to the original, deviance-based selection suggested for componentwise boosting. With \code{score}, a score statistic is employed. The latter is faster and corresponds to the approach employed for Cox regression in package \code{CoxBoost}.}
\item{stepsize.factor.linear}{determines the step-size modification factor by which the natural step size of boosting steps for the linear covariates should be changed after a covariate has been selected in a boosting step. The default (value \code{1}) implies constant penalties, for a value < 1 the penalty for a covariate is increased after it has been selected in a boosting step, and for a value > 1 the penalty it is decreased. If \code{pendistmat} is given, penalty updates are only performed for covariates that have at least one connection to another covariate.}
\item{sf.scheme}{scheme for changing step sizes (via \code{stepsize.factor}). \code{"linear"} corresponds to the scheme described in Binder and Schumacher (2009), \code{"sigmoid"} employs a sigmoid shape.}
\item{pendistmat.linear}{connection matrix with entries ranging between 0 and 1, with entry \code{(i,j)} indicating the certainty of the connection between covariates \code{i} and \code{j}. According to this information penalty changes due to \code{stepsize.factor.linear} < 1 are propagated, i.e., if entry \code{(i,j)} is non-zero, the penalty for covariate \code{j} is decreased after it has been increased for covariate \code{i}, after it has been selected in a boosting step. This matrix either has to have dimension \code{q * q} or the indicices of the \code{q.connected} connected linear covariates have to be given in \code{connected.index.linear}, in which case the matrix has to have dimension \code{q.connected * q.connected}.}
\item{connected.index.linear}{indices of the \code{q.connected} connected linear covariates, for which \code{pendistmat.linear} provides the connection information for distributing changes in penalties. If \code{NULL}, and a connection matrix is given, all covariates are assumed to be connected.}
\item{weights}{an optional vector of weights to be used in the fitting process.}
\item{stepno}{number of boosting steps (\code{m}).}
\item{family}{a description of the error distribution to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) Note that \code{GAMBoost} supports only canonical link functions and no scale parameter estimation, so \code{gaussian()}, \code{binomial()} and \code{poisson()} are the most plausible candidates.}
\item{sparse.boost}{logical value indicating whether a criterion considering degrees of freedom (specifically AIC) should be used for selecting a covariate for an update in each boosting step (sparse boosting), instead of the deviance.}
\item{sparse.weight}{factor modifying how the degrees of freedom enter into calculation of AIC for sparse boosting.}
\item{calc.hat}{logical value indicating whether the hat matrix should be computed for each boosting step. If set to \code{FALSE} no degrees of freedom and therefore e.g. no AIC will be available. On the other hand fitting will be faster (especially for a large number of observations).}
\item{calc.se}{logical value indicating whether confidence bands should be calculated. Switch this of for faster fitting.}
\item{AIC.type}{type of model selection criterion to be calculated (and also to be used for sparse boosting if applicable): in the Gaussian case \code{"classical"} gives AIC and \code{"corrected"} results in corrected AIC (Hurvich, Simonoff and Tsai, 1998); for all other response families standard AIC is used.}
\item{return.score}{logical value indicating whether the value of the score statistic, as evaluated in each boosting step for every covariate (only for binary response models, for linear components, if \code{criterion="score"}), should be returned. The corresponding element \code{scoremat} can become very large (and needs much memory) when the number of covariates and boosting steps is large.}
\item{trace}{logical value indicating whether progress in estimation should be indicated by printing the index  of the covariate updated (\code{1}, \dots , \code{p} for smooth components and \code{p+1}, \dots ,\code{p-q} for parametric components) in the current boosting step.}
}
\details{
The idea of likelihood based boosting (Tutz and Binder, 2006) is most easily understood for models with a Gaussian response.  There it results in repeated fitting of residuals (This idea is then transferred to the generalized case).  For obtaining an additive model GAMBoost uses a large number of boosting steps where in each step a penalized B-spline (of degree \code{bdeg}) (Eilers and Marx, 1996) is fitted to one covariate, the response being the residuals from the last step.  The covariate to be updated is selected by deviance (or in case of sparse boosting by some model selection criterion).  The B-spline coefficient estimates in each step are fitted under the constraint of a large penalty on their first (or higher order) differences. So in each step only a small adjustment is made.  Summing over all steps for each covariate a smooth function estimate is obtained. When no basis expansion is used, i.e. just coefficients of covariates are updated, (generalized) linear models are obtained.  

The main parameter of the algorithm is the number of boosting steps, given that the penalty is chosen large enough (If too small the minimum AIC will occur for a very early boosting step; see \code{\link{optimGAMBoostPenalty}}).  When there is a large number of covariates with potentially non-linear effect, having a single parameter (with adaptive smoothness assignment to single components performed automatically by the algorithm) is a huge advantage compared to approaches, where a smoothing parameter has to be selected for each single component. The biggest advantage over conventional methods for fitting generalized additive models (e.g. \code{mgcv:gam} or \code{gam:gam}) will therefore be obtained for a large number of covariates compared to the number of observations (e.g. 10 covariates with only 100 observations).  In addition GAMBoost performs well compared to other approaches when there is a small signal-to-noise ratio and/or a response (e.g. binary) with a small amount of information.

If a group of correlated covariates has influence on the response, e.g. genes from the same pathway, componentwise boosting will often result in a non-zero estimate for only one member of this group. To avoid this, information on the connection between covariates with linear influence can be provided in \code{pendistmat.linear}. If then, in addition, a penalty updating scheme with \code{stepsize.factor.linear} < 1 is chosen, connected covariates with linear influence are more likely to be chosen in future boosting steps, if a directly connected covariate has been chosen in an earlier boosting step (see Binder and Schumacher, 2008b).

Note that the degrees of freedom (and based on these AIC and BIC) are just approximations, which are completely valid for example only when the order and the indices of the components updated is fixed.  This leads to problems especially when there is a very large number of covariates (e.g. 10000 covariates with only 100 observations).  Then it might be better (but also slower) to rely on cross validation (see \code{\link{cv.GAMBoost}}) for selection of the number of boosting steps.  

Note that the \code{gamboost} routine in the R package \code{mboost} implements a different kind of boosting strategy: gradient based boosting instead of likelihood based boosting.  The two approaches coincide only in special cases (e.g. L2 loss  and Gaussian response).  While gradient based boosting is more general, only likelihood based boosting allows e.g. for easily obtainable pointwise confidence bands.  
}
\value{
\code{GAMBoost} returns an object of class \code{GAMBoost}.  \code{GAMBoost} objects can be examined by print, summary, and plot.  \code{\link{getGAMBoostSelected}} can be used on them to identify selected/significant covariates.  

\item{x, n, p.linear}{original values for covariates with non-linear effect, number of observations, and number of covariates with linear effect.}
\item{penalty, penalty.linear}{penalties used in updating smooth and linear components.}
\item{stepno}{number of boosting steps.}
\item{family}{response family.}
\item{AIC.type}{type of AIC given in component \code{AIC} (applies only for Gaussian response).}
\item{deviance, trace, AIC, BIC}{vectors of length \code{m} giving deviance, approximate degrees of freedom, AIC, and BIC for each boosting step.}
\item{selected}{vector of length \code{m} given the index of the covariate updated in each boosting step (\code{1-p} for smooth components and \code{(p+1)}, ... , \code{(p-q)} for parametric components).}
\item{beta}{list of length \code{p+1}, each element being a matrix with \code{m+1} rows giving the estimated coefficients for the intercept term (\code{beta[[1]]}) and the smooth terms (\code{beta[[2]]}, ... , \code{beta[[p+1]]}).}
\item{beta.linear}{\code{m * q} matrix containing the coefficient estimates for the (standardized) linear covariates.}
\item{scoremat}{\code{m * q} matrix containing the value of the score statistic for each of the optional covariates before each boosting step, if selection has been performed by the score statistic (\code{criterion="score"}) and a binary response model has been employed.}
\item{mean.linear, sd.linear}{vector of mean values and standard deviations used for standardizing the linear covariates.}
\item{hatmatrix}{hat matrix at the final boosting step.}
\item{eta}{\code{n * (m+1)} matrix with predicted value (at predictor level) for each boosting step.}
\item{predictors}{list of length \code{p+1} containing information (as a list) on basis expansions for the smooth model components (\code{predictors[[2]]}, ... , \code{predictors[[p+1]]}).}
}
\author{
Written by Harald Binder \email{binderh@fdm.uni-freiburg.de}, matching closely the original Fortran implementation employed
for Tutz and Binder (2006). 
}
\references{
Binder, H. and Schumacher, M. (2009). Incorporating pathway information into boosting estimation of high-dimensional risk prediction models. BMC Bioinformatics. 10:18.

Binder, H. and Schumacher, M. (2008). Incorporating pathway information into boosting estimation of high-dimensional risk prediction models. Manuscript.

Hurvich, C. M., Simonoff, J. S. and Tsai, C. L. (1998). Smoothing parameter selection in nonparametric regression using and improved Akaike information criterion. Journal of the Royal Statistical Society B, \bold{60}(2), 271--293.

Eilers, P. H. C. and Marx, B. D. (1996) Flexible smoothing with B-splines and penalties. \emph{Statistical Science}, \bold{11}(2), 89--121.

Tutz, G. and Binder, H. (2007) Boosting ridge regression. \emph{Computational Statistics \& Data Analysis}, \bold{51}(12), 6044--6059.

Tutz, G. and Binder, H. (2006) Generalized additive modelling with implicit variable selection by likelihood based boosting. \emph{Biometrics}, \bold{51}, 961--971.
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
plot(0:100,gb3$AIC,type="l",xlab="stepno",ylab="AIC"); lines(0:100,gb4$AIC,lty=2)


##   Fit a generalized linear model with penalty modification
##   after every boosting step, with penalty changes being
##   redistrbuted via a connection matrix

pendistmat <- matrix(0,10,10)
#    Covariates 1 and 3 are connected
pendistmat[1,3] <- pendistmat[3,1] <- 1

gb5 <- GAMBoost(x=NULL,y=y,x.linear=x.linear,penalty.linear=100,
                stepsize.factor.linear=0.9,pendistmat.linear=pendistmat,
                stepno=100,trace=TRUE,family=binomial())

#   or alternatively

gb5 <- GAMBoost(x=NULL,y=y,x.linear=x.linear,penalty.linear=100,
                stepsize.factor.linear=0.9,
                pendistmat.linear=pendistmat[c(1,3),c(1,3)],
                connected.index.linear=c(1,3),
                stepno=100,trace=TRUE,family=binomial())

}
\keyword{models} \keyword{smooth} \keyword{regression}
