GAMBoost <- function(x=NULL,y,xmin=NULL,xmax=NULL,penalty=100,bdeg=2,pdiff=1,
                     x.linear=NULL,standardize.linear=TRUE,penalty.linear=0,
                     weights=rep(1,length(y)),stepno=500,family=binomial(),
                     sparse.boost=FALSE,sparse.weight=1,calc.hat=TRUE,calc.se=TRUE,
                     AIC.type="corrected",trace=FALSE) 
{   
    #
    #   check parameters for consistency and transform them to the format expected
    #
     
    if (sparse.boost) calc.hat <- TRUE
    if (calc.se) calc.hat <- TRUE
    
    #   get argument 'family' into the format needed
    if (is.character(family)) {
        family <- switch(family,gaussian=gaussian(),binomial=binomial(),poisson=poisson())
    } else {
        if (is.function(family)) family <- family()
    }
    
    if (sum(weights) != length(y)) {
        weights <- weights / (sum(weights)/length(y))   
    }
      
    #   consistency checks
    
    if (is.null(x) && is.null(x.linear)) {
        cat("Neither non-parametric components (x) nor parametric components (x.linear) given. Abort.\n")
        return(NULL)
    }
    
    canonical.link <- switch(family$family,gaussian="identity",binomial="logit",poisson="log") 
    if (is.null(canonical.link) || family$link != canonical.link) {
        warning(paste("GAMBoost expects the canonical link for family '",family$family,"' and does not honor dispersion parameters.\n",sep=""))
    }
    
    #
    #   initialized data structures
    #
    
    predictors <- basis.expansion(x,n=length(y),xmin=xmin,xmax=xmax,bdeg=bdeg,pdiff=pdiff)
    
    #   If fstruct gets too big some time then inline this function call
    #   For now it makes the code cleaner
    fstruct <- init.fit.structure(predictors,stepno=stepno,family=family,calc.se=calc.se)
    fstruct$penalty <- penalty
    fstruct$AIC.type <- AIC.type
    fstruct$x <- x


    #   initialize/standardize linear predictors (if given)

    if (!is.null(x.linear)) {
        if (standardize.linear) {
            fstruct$mean.linear <- apply(x.linear,2,mean)
            fstruct$sd.linear <- apply(x.linear,2,sd)
        } else {
            fstruct$mean.linear <- rep(0,ncol(x.linear))
            fstruct$sd.linear <- rep(1,ncol(x.linear))
        }
        
        fstruct$x.linear <- t((t(x.linear) - fstruct$mean.linear)/fstruct$sd.linear)
        fstruct$beta.linear <- matrix(0,fstruct$stepno+1,ncol(fstruct$x.linear))   
        
        if (length(penalty.linear) < ncol(x.linear)) penalty.linear <- rep(penalty.linear[1],ncol(x.linear))
    } else {
        fstruct$x.linear <- NULL
        fstruct$beta.linear <- NULL
    }
    
    fstruct$penalty.linear <- penalty.linear
    
    
    #   
    #   Boosting iterations
    #
    
    for (actual.step in 1:stepno) {
        #cat("*** iteration",actual.step,"\n")
        
        #   update intercept
        
        intercept.D <- fstruct$family$mu.eta(fstruct$eta[,actual.step])
        
        if (!is.null(fstruct$x.linear) && sum(fstruct$penalty.linear == 0) > 0) {
            #   update unpenalized linear covariates together with intercept 
            
            unpen.x <- cbind(rep(1,fstruct$n),fstruct$x.linear[,fstruct$penalty.linear == 0])
            unpen.pre <- solve(t((unpen.x*intercept.D)*weights) %*% unpen.x) %*% t(unpen.x)
                                     
            unpen.beta.delta <- as.vector(unpen.pre %*% (weights*(y - fstruct$family$linkinv(fstruct$eta[,actual.step]))))
            unpen.eta.delta <- unpen.x %*% unpen.beta.delta

            fstruct$beta[[1]][actual.step+1,] <- fstruct$beta[[1]][actual.step,] + unpen.beta.delta[1]
            fstruct$beta.linear[actual.step+1,fstruct$penalty.linear == 0] <- fstruct$beta.linear[actual.step,fstruct$penalty.linear == 0] + unpen.beta.delta[2:length(unpen.beta.delta)]
            fstruct$eta[,actual.step+1] <- fstruct$eta[,actual.step] + unpen.x %*% unpen.beta.delta
            
            if (calc.hat) {
                pre.hat <- (unpen.x %*% unpen.pre) * intercept.D             
                if (actual.step == 1) {
                    fstruct$hatmatrix[] <- pre.hat
                } else {
                    fstruct$hatmatrix <- fstruct$hatmatrix + pre.hat - pre.hat %*% fstruct$hatmatrix
                }
            }
        } else {
            #   update just the untercept
            
            intercept.pre <- 1/sum(intercept.D*weights)
            intercept.beta.delta <- sum(intercept.pre * (weights*(y - fstruct$family$linkinv(fstruct$eta[,actual.step]))))
            fstruct$beta[[1]][actual.step+1,] <- fstruct$beta[[1]][actual.step,] + intercept.beta.delta
            fstruct$eta[,actual.step+1] <- fstruct$eta[,actual.step] + intercept.beta.delta

            if (calc.hat) {
                pre.hat <- intercept.pre*intercept.D
                if (actual.step == 1) {
                    fstruct$hatmatrix[] <- pre.hat
                } else {
                    #   efficient update using the special structure of the intercept pre.hat
                    fstruct$hatmatrix <- fstruct$hatmatrix + matrix(1 - apply(fstruct$hatmatrix,2,sum),fstruct$n,fstruct$n,byrow=TRUE) * pre.hat
                }
            }
        }

        
        #   update estimates for penalized covariates

        actual.eta <- fstruct$eta[,actual.step+1]
        actual.mu <- as.vector(fstruct$family$linkinv(actual.eta))
        D <- fstruct$family$mu.eta(actual.eta)

        #   the following works only if the canonical link is used

        best.candidate <- 1
        best.deviance <- -1
        best.criterion <- -1
        best.pre <- NULL
        best.beta.delta <- NULL
        best.eta <- NULL
       
        if (sparse.boost) actual.trace <- sum(diag(fstruct$hatmatrix))
       
        #   
        #   evaluate non-parametric components
        #
        
        if (length(predictors) > 1) {
            for (actual.predictor in 2:length(predictors)) {
                #cat("evaluating",actual.predictor,"\n")

                candidate.pre <- solve(t((predictors[[actual.predictor]]$expansion*D)*weights) %*% 
                                         predictors[[actual.predictor]]$expansion + 
                                         penalty*predictors[[actual.predictor]]$penalty) %*% 
                                         t(predictors[[actual.predictor]]$expansion)
                                 
                candidate.beta.delta <- as.vector(candidate.pre %*% (weights*(y - actual.mu)))
                candidate.eta.delta <- predictors[[actual.predictor]]$expansion %*% candidate.beta.delta

                candidate.mu <- fstruct$family$linkinv(fstruct$eta[,actual.step+1] + candidate.eta.delta)
                candidate.deviance <- sum(fstruct$family$dev.resids(y,candidate.mu,weights))
        
                #   sparse boosting used AIC instead of deviance to select a predictor for update
                if (sparse.boost) {
                    candidate.pre.hat <- (predictors[[actual.predictor]]$expansion %*% candidate.pre) * D
                    #   we don't need the whole candidate hat matrix, just the trace
                    candidate.trace <- actual.trace + sum(diag(candidate.pre.hat)) - sum(candidate.pre.hat * t(fstruct$hatmatrix))
                    if (fstruct$family$family == "gaussian") {
                        if (fstruct$AIC.type == "corrected") {
                            candidate.criterion <- log(candidate.deviance/fstruct$n) + (1+sparse.weight*candidate.trace/fstruct$n)/(1-(sparse.weight*candidate.trace+2)/fstruct$n)  
                        } else {
                            candidate.criterion <- log(candidate.deviance/fstruct$n) + sparse.weight*2*candidate.trace/fstruct$n  
                        }
                    } else {
                        candidate.criterion <- candidate.deviance + sparse.weight*2*candidate.trace  
                    }
                } else {
                    candidate.criterion <- candidate.deviance
                }      

                if (actual.predictor == 2 || candidate.criterion < best.criterion) {
                    best.candidate <- actual.predictor
                    best.deviance <- candidate.deviance
                    best.criterion <- candidate.criterion
                    best.pre <- candidate.pre
                    best.beta.delta <- candidate.beta.delta
                    best.eta <- fstruct$eta[,actual.step+1] + candidate.eta.delta
                }
            }
        }
        
        #
        #   evaluate penalized linear components
        #
        
        if (!is.null(fstruct$x.linear) && sum(fstruct$penalty.linear != 0) > 0) {
            
            pre.mult <- (1/(apply(fstruct$x.linear[,fstruct$penalty.linear != 0,drop=FALSE]^2 * D * weights,2,sum) + fstruct$penalty.linear[fstruct$penalty.linear != 0]))  
            pre.mat <- t(fstruct$x.linear[,fstruct$penalty.linear != 0,drop=FALSE]) * pre.mult
                            
              
            beta.delta.vec <- drop(pre.mat %*% (weights*(y - actual.mu)))
 
            dev.vec <- apply(matrix(
                                family$dev.resids(rep(y,sum(fstruct$penalty.linear != 0)),
                                    family$linkinv(as.vector(t(t(fstruct$x.linear[,fstruct$penalty.linear != 0,drop=FALSE]) * beta.delta.vec) + actual.eta))
                                    ,rep(weights,sum(fstruct$penalty.linear != 0))),
                                nrow(fstruct$x.linear),sum(fstruct$penalty.linear != 0)),2,sum)
            
            crit.vec <- dev.vec

            if (sparse.boost) {
                pen.lin.index <- (1:length(fstruct$penalty.linear))[fstruct$penalty.linear != 0]
                
                for (i in 1:length(pen.lin.index)) {
                    candidate.pre.hat <- (fstruct$x.linear[,pen.lin.index[i],drop=FALSE] %*% pre.mat[i,,drop=FALSE]) * D
                    #   we don't need the whole candidate hat matrix, just the trace
                    candidate.trace <- actual.trace + sum(fstruct$x.linear[,pen.lin.index[i]]^2*D)*pre.mult[i] - sum(candidate.pre.hat * t(fstruct$hatmatrix))
                    
                    if (fstruct$family$family == "gaussian") {
                        if (fstruct$AIC.type == "corrected") {
                            crit.vec[i] <- log(crit.vec[i]/fstruct$n) + (1+sparse.weight*candidate.trace/fstruct$n)/(1-(sparse.weight*candidate.trace+2)/fstruct$n)  
                        } else {
                            crit.vec[i] <- log(crit.vec[i]/fstruct$n) + sparse.weight*2*candidate.trace/fstruct$n  
                        }
                    } else {
                        crit.vec[i] <- crit.vec[i] + sparse.weight*2*candidate.trace  
                    }
                }
            }
                                
            best.candidate.linear <- (1:length(fstruct$penalty.linear))[fstruct$penalty.linear != 0][which.min(crit.vec)]
            
            if (best.criterion == -1 || min(crit.vec) < best.criterion) {
                best.candidate <- best.candidate.linear + length(predictors)
                best.deviance <- min(dev.vec)
                best.criterion <- min(crit.vec)
                best.pre <- 1/((fstruct$x.linear[,best.candidate.linear]*D*weights) %*% fstruct$x.linear[,best.candidate.linear] + fstruct$penalty.linear[best.candidate.linear]) * fstruct$x.linear[,best.candidate.linear]
                best.beta.delta <- best.pre %*% (weights*(y - actual.mu))
                best.eta <- actual.eta + fstruct$x.linear[,best.candidate.linear] * best.beta.delta
            }
        }

        #cat("updating",best.candidate-1,"with deviance",best.deviance,"\n")
        if (trace) cat(best.candidate-1," ")
        
        fstruct$selected <- c(fstruct$selected,best.candidate-1)
        fstruct$eta[,actual.step+1] <- as.vector(best.eta)
        
        if (length(predictors) > 1) {
            for (j in 2:length(predictors)) {
                if (j == best.candidate) {
                    fstruct$beta[[j]][actual.step+1,] <- fstruct$beta[[j]][actual.step,] + as.vector(best.beta.delta)
                } else {
                    fstruct$beta[[j]][actual.step+1,] <- fstruct$beta[[j]][actual.step,]
                }
            }
        }
        
        if (!is.null(fstruct$x.linear) && sum(fstruct$penalty.linear != 0) > 0) {
            fstruct$beta.linear[actual.step+1,fstruct$penalty.linear != 0] <- fstruct$beta.linear[actual.step,fstruct$penalty.linear != 0]
            if (best.candidate > length(predictors)) fstruct$beta.linear[actual.step+1,best.candidate - length(predictors)] <- fstruct$beta.linear[actual.step+1,best.candidate - length(predictors)] + best.beta.delta
        }
        
        if (calc.hat || calc.se) {
            if (best.candidate <= length(predictors)) {
                pre.Q <- (predictors[[best.candidate]]$expansion %*% best.pre)
    
                if (calc.se) {
                    #   check whether the confidence bands for this covariate have to be initialized
                    if (is.na(fstruct$Qmatrix[[best.candidate]][1])) {
                        fstruct$Qmatrix[[best.candidate]] <- matrix(0,fstruct$n,fstruct$n)
                        fstruct$obsvar[[best.candidate]] <- matrix(NA,fstruct$n,stepno+1)
                    }
                
                    fstruct$Qmatrix[[best.candidate]] <- fstruct$Qmatrix[[best.candidate]] + pre.Q - pre.Q %*% fstruct$hatmatrix
                    fstruct$obsvar[[best.candidate]][,actual.step + 1] <- diag(fstruct$Qmatrix[[best.candidate]] %*% 
                                                                              (t(fstruct$Qmatrix[[best.candidate]]) * fstruct$family$variance(fstruct$family$linkinv(fstruct$eta[,actual.step+1]))))
                }
            } else {
                pre.Q <- (fstruct$x.linear[,best.candidate-length(predictors),drop=FALSE] %*% best.pre)
            }
    
            fstruct$hatmatrix <- fstruct$hatmatrix + (pre.Q * D) - (pre.Q * D) %*% fstruct$hatmatrix
            fstruct$trace <- c(fstruct$trace,sum(diag(fstruct$hatmatrix)))
        }
        
        fstruct$deviance <- c(fstruct$deviance,best.deviance)
    }
    if (trace) cat("\n")
    
    if (calc.hat) {
        if (fstruct$family$family == "gaussian") {
            #   estimate sigma^2 from the largest model
            sigma.sq.hat <- fstruct$deviance[length(fstruct$deviance)]/(fstruct$n-fstruct$trace[length(fstruct$trace)])
            if (fstruct$AIC.type == "corrected") {
                fstruct$AIC <- log(fstruct$deviance/fstruct$n) + (1+fstruct$trace/fstruct$n)/(1-(fstruct$trace+2)/fstruct$n)
            } else {
                fstruct$AIC <- fstruct$deviance/sigma.sq.hat + 2*(fstruct$trace+1) 
            }
            fstruct$BIC <- fstruct$deviance/sigma.sq.hat + log(fstruct$n)*(fstruct$trace+1)
        } else {
            fstruct$AIC <- fstruct$deviance + 2*fstruct$trace
            fstruct$BIC <- fstruct$deviance + log(fstruct$n)*fstruct$trace
        }
    }
    
    fstruct$predictors <- predictors
    
    #   remove some intermediate structures for saving space
    
    fstruct$Qmatrix <- NULL
    
    class(fstruct) <- "GAMBoost"
    return(fstruct)
}

basis.expansion <- function(x,n,xmin=NULL,xmax=NULL,bdeg=2,ndx=20,pdiff=1) {
    if (is.null(xmin) && !is.null(x)) xmin <- apply(x,2,min)
    if (is.null(xmax) && !is.null(x)) xmax <- apply(x,2,max)

    result <- list()
    
    #   intercept term
    result[[1]] <- list(expansion=matrix(rep(1,n),n,1))
    
    #   currently all predictors get a B-spline expansion    
    if (!is.null(x)) {
        for (i in 1:ncol(x)) {
            predictor <- list()
            predictor$xl <- xmin[i]
            predictor$xr <- xmax[i]
            predictor$bdeg <- bdeg
            dx <- (predictor$xr - predictor$xl)/ndx
            predictor$knots <- seq(predictor$xl - predictor$bdeg * dx, predictor$xr + predictor$bdeg * dx, by = dx)
            predictor$knots[bdeg+1] <- predictor$xl
            predictor$knots[length(predictor$knots) - predictor$bdeg] <- predictor$xr
            predictor$expansion <- spline.des(predictor$knots, x[,i], predictor$bdeg + 1, rep(0,nrow(x)))$design
            predictor$expansion <- predictor$expansion %*% rbind(diag(ncol(predictor$expansion)-1),rep(-1,ncol(predictor$expansion)-1))
        
            k <- diff(rbind(diag(ncol(predictor$expansion)),rep(-1,ncol(predictor$expansion))),,pdiff)
            predictor$penalty <- t(k)%*%k
        
            result[[length(result)+1]] <- predictor
        }
    }

    return(result)
}

init.fit.structure <- function(predictors,stepno,family,calc.se) {
    fit.structure <- list()
    fit.structure$n <- nrow(predictors[[1]]$expansion)
    fit.structure$stepno <- stepno
    fit.structure$family <- family
    fit.structure$eta <- matrix(0,fit.structure$n,stepno+1)
    fit.structure$hatmatrix <- matrix(0,fit.structure$n,fit.structure$n)
    fit.structure$selected <- c()
    fit.structure$deviance <- c()
    fit.structure$trace <- c()
    
    beta <- list()
    beta[[1]] <- matrix(0,stepno+1,1)
    Qmatrix <- list()
    Qmatrix[[1]] <- NA
    obsvar <- list()
    obsvar[[1]] <- NA
    
    if (length(predictors) > 1) {
        for (i in 2:length(predictors)) {
            beta[[i]] <- matrix(0,stepno+1,ncol(predictors[[i]]$expansion))
            if (calc.se) {
                Qmatrix[[i]] <- NA
                obsvar[[i]] <- NA

                #Qmatrix[[i]] <- matrix(0,fit.structure$n,fit.structure$n)
                #obsvar[[i]] <- matrix(NA,fit.structure$n,stepno+1)
            } else {
                Qmatrix[[i]] <- NA
                obsvar[[i]] <- NA
            }
        }
    }
    
    fit.structure$beta <- beta
    fit.structure$Qmatrix <- Qmatrix
    fit.structure$obsvar <- obsvar
    
    return(fit.structure)
}

getGAMBoostSelected <- function(object,at.step=NULL) {
    if (is.null(at.step)) at.step <- object$stepno
    
    result <- list(smooth=c(),smoothbands=c(),parametric=c())
    
    if (length(object$predictors) > 1) {
        for (i in 2:length(object$beta)) {
            if (sum(object$beta[[i]][at.step+1,] != 0) > 0) result$smooth <- c(result$smooth,i-1)            
        }
        if (length(object$obsvar[[1+1]]) > 1) { #   variance information available
            for (i in 2:length(object$beta)) {
                band.info <- null.in.bands(object,select=i-1,at.step=at.step)
                if (!is.na(band.info) && !band.info) result$smoothbands <- c(result$smoothbands,i-1)
            }
        }
    }

    if (!is.null(object$x.linear)) {
         for (i in 1:ncol(object$x.linear)) {
             if (object$beta.linear[at.step+1,i] != 0) result$parametric <- c(result$parametric,i)
         }
    }

    return(result)
}

summary.GAMBoost <- function(object,...) {
    print(object)
    
    cat("\nfitted covariates:\n")

    fitted.covariates <- function(at.step) {
        
        retstring <- paste("(step ",at.step,"):\n",sep="")
        
        retstring <- paste(retstring,"        intercept: ",round(object$beta[[1]][at.step+1,1],4),"\n",sep="")
        
        selected <- getGAMBoostSelected(object,at.step=at.step)
        
        if (length(object$predictors) > 1) {
            retstring <- paste(retstring,"        smooth selected: ",sep="")
        
            smooth.names <- colnames(object$x)
            if (is.null(smooth.names)) smooth.names <- paste("S",1:(length(object$predictors)-1),sep="")
            
            if (length(selected$smooth) > 0) retstring <- paste(retstring,paste(smooth.names[selected$smooth],collapse=", "))
            retstring <- paste(retstring,"\n",sep="")
            
            if (length(object$obsvar[[1+1]]) > 1) { #   variance information available
                retstring <- paste(retstring,"        smooth bands don't contain zero: ",sep="")
                if (length(selected$smoothbands) > 0) retstring <- paste(retstring,paste(smooth.names[selected$smoothbands],collapse=", "))
                retstring <- paste(retstring,"\n",sep="")
            }
        }

        if (!is.null(object$x.linear)) {
            retstring <- paste(retstring,"        parametric selected: ",sep="")
            
            parametric.names <- colnames(object$x.linear)
            if (is.null(parametric.names)) parametric.names <- paste("V",1:ncol(object$x.linear),sep="")
            if (length(selected$parametric) > 0) retstring <- paste(retstring,paste(parametric.names[selected$parametric]," (",round(object$beta.linear[at.step+1,selected$parametric],4),")",collapse=", ",sep=""),sep="")
            retstring <- paste(retstring,"\n",sep="")
            
       }
    
        return(retstring)
    }

    cat("    at final boosting step",fitted.covariates(object$stepno),"\n")
    if (!is.null(object$trace)) {
        cat("    at minimum AIC ",ifelse(object$family$family=="gaussian",paste("(",object$AIC.type,") ",sep=""),""),fitted.covariates(which.min(object$AIC)),"\n",sep="")        
        cat("    at minimum BIC",fitted.covariates(which.min(object$BIC)),"\n")        
    }
}

print.GAMBoost <- function(x,...) {
    cat("\n")
    cat("family:",x$family$family,"(with canonical link)\n")
    cat("model components: ",length(x$predictors)-1," smooth",
        ifelse(length(x$predictors) > 1,paste(" (with penalty ",x$penalty,")",sep=""),""),", ",
        ifelse(!is.null(x$x.linear),ncol(x$x.linear),0)," parametric",
        ifelse(!is.null(x$x.linear),paste(" (with penalty ",
            ifelse(sum(x$penalty.linear == mean(x$penalty.linear)) != length(x$penalty.linear),
                paste(x$penalty.linear,collapse=", "),x$penalty.linear[1]),")",sep=""),""),"\n",sep="")
    
    cat("\nmodel fit:\n")
    
    mfit.string <- function(at.step) return(paste("(step ",at.step,"):\n        residual deviance ",
                                            round(x$deviance[at.step],4),ifelse(!is.null(x$trace),
                                            paste(", df ",round(x$trace[at.step],4),", AIC ",ifelse(x$family$family=="gaussian",paste("(",x$AIC.type,") ",sep=""),""),round(x$AIC[at.step],4),", BIC ",round(x$BIC[at.step],4),sep=""),""),sep=""))
    
    cat("    at final boosting step",mfit.string(x$stepno),"\n")
    if (!is.null(x$trace)) {
        cat("    at minimum AIC ",ifelse(x$family$family=="gaussian",paste("(",x$AIC.type,") ",sep=""),""),mfit.string(which.min(x$AIC)),"\n",sep="")        
        cat("    at minimum BIC",mfit.string(which.min(x$BIC)),"\n")        
    }
}

predict.GAMBoost <- function(object,newdata=NULL,newdata.linear=NULL,at.step=NULL,type="link",...) {
    if (is.null(at.step)) { 
        at.step <- nrow(object$beta[[1]])
    } else {
        at.step <- at.step + 1
    }
    
    n <- ifelse(is.null(newdata),object$n,nrow(newdata))
    
    if (is.null(newdata) && type != "terms") {
        eta <- object$eta[,at.step]
    } else {
        eta <- matrix(rep(object$beta[[1]][at.step,],n),n,length(at.step),byrow=TRUE)
        
        if (length(object$predictors) > 1) {
            for (i in 2:length(object$predictors)) {
                if (is.null(newdata)) {
                    actual.expansion <- object$predictors[[i]]$expansion
                } else {
                    actual.expansion <- spline.des(object$predictors[[i]]$knots, newdata[,i-1], 
                                                   object$predictors[[i]]$bdeg + 1, rep(0,nrow(newdata)),outer.ok=TRUE)$design
                    actual.expansion <- actual.expansion %*% rbind(diag(ncol(actual.expansion)-1),rep(-1,ncol(actual.expansion)-1))
                }
            
                if (type == "terms") {
                    eta <- cbind(eta,actual.expansion %*% object$beta[[i]][at.step,])                
                } else {
                    eta <- eta + actual.expansion %*% t(object$beta[[i]][at.step,,drop=FALSE])
                }
            }
        }

        if (type != "terms" && !is.null(object$x.linear)) {
            if (!is.null(newdata.linear)) {
                eta <- eta + t((t(newdata.linear) - object$mean.linear) / object$sd.linear) %*% t(object$beta.linear[at.step,,drop=FALSE])
            } else {
                warning("No data for parametric components given. Their contribution for prediction is ignored.\n")
            }
        }
    }
    
    if (length(eta) == n) eta <- drop(eta)
    
    if (type == "terms") return(eta[,2:ncol(eta)])    
    if (type == "link") return(eta)
    
    return(object$family$linkinv(eta))
}

plot.GAMBoost <- function(x,select=NULL,at.step=NULL,add=FALSE,phi=1,ylim=NULL,xlab=NULL,ylab=NULL,...) {
    if (length(x$predictors) < 2) {
        cat("No non-paramatric components fitted, i.e. no plots available.\n")
        return(NULL)
    }
    if (is.null(select)) select <- 1:(length(x$predictors)-1)
    if (is.null(at.step)) { 
        at.step <- nrow(x$beta[[1]])
    } else {
        at.step <- at.step + 1
    }

    smooth.names <- colnames(x$x)
    if (is.null(smooth.names)) smooth.names <- paste("S",1:(length(x$predictors)-1),sep="")

    set.xlab <- is.null(xlab)

    for (i in select) {
        actual.x <- seq(from=x$predictors[[i+1]]$xl,to=x$predictors[[i+1]]$xr,length=100)
        actual.expansion <- spline.des(x$predictors[[i+1]]$knots, actual.x, x$predictors[[i+1]]$bdeg + 1, rep(0,length(actual.x)))$design
        actual.expansion <- actual.expansion %*% rbind(diag(ncol(actual.expansion)-1),rep(-1,ncol(actual.expansion)-1))
        eta <- as.vector(actual.expansion %*% x$beta[[i+1]][at.step,])
        
        bands <- calc.confidence.bands(x,i,at.step-1,phi)
        
        if (!is.null(bands)) {
            if (is.null(ylim)) ylim <- c(min(bands$lower),max(bands$upper))            
        }

        if (is.null(ylim)) ylim <- c(min(eta),max(eta))
        
        if (!add) {
            if (set.xlab) xlab <- smooth.names[i]
            if (is.null(ylab)) ylab <- "eta"
            plot(actual.x,eta,type="l",ylim=ylim,xlab=xlab,ylab=ylab,...)
        } else {
            lines(actual.x,eta)
        }

        if (!is.null(bands)) {
            lines(bands$x,bands$upper,lty=2)
            lines(bands$x,bands$lower,lty=2)            
        }
    }
}

calc.confidence.bands <- function(object,covno,at.step,phi) {
    at.step <- at.step + 1
    result <- NULL
    if (length(object$obsvar[[covno+1]]) > 1) {    #   there is variance information available
        updated <- (1:at.step)[!is.na(object$obsvar[[covno+1]][1,1:at.step])]
        if (length(updated) > 0) {
            actual.obsvar <- object$obsvar[[covno+1]][,max(updated)]
            ori.x <- object$x[,covno]
            ori.eta <- predict(object,type="terms",at.step=at.step-1)[,covno]

            unique.x <- sort(unique(ori.x))
            mean.eta <- rep(0,length(unique.x))
            upper.eta <- rep(0,length(unique.x))
            lower.eta <- rep(0,length(unique.x))

            for (j in 1:length(unique.x)) {
                this.sd <- sqrt(mean(actual.obsvar[ori.x == unique.x[j]],na.rm=TRUE))
                mean.eta[j] <- mean(ori.eta[ori.x == unique.x[j]])
                upper.eta[j] <- mean.eta[j] + 2*this.sd * sqrt(phi)
                lower.eta[j] <- mean.eta[j] - 2*this.sd * sqrt(phi)
            }
     
            result <- list(x=unique.x,eta=mean.eta,upper=upper.eta,lower=lower.eta)
        }
    }
    return(result)    
}

null.in.bands <- function(object,select=NULL,at.step=NULL,phi=1) {
    if (is.null(select)) select <- 1:(length(object$predictors)-1)
    if (is.null(at.step)) { 
        at.step <- nrow(object$beta[[1]]) - 1
    } 

    result <- matrix(NA,length(at.step),length(select))
    for (i in 1:length(select)) {
        for (j in 1:length(at.step)) {
            bands <- calc.confidence.bands(object,select[i],at.step[j],phi)
            if (!is.null(bands)) {
                if (sum(bands$upper < 0) + sum(bands$lower > 0) == 0) {
                    result[j,i] <- TRUE
                } else {
                    result[j,i] <- FALSE
                }
            }
        }
    }
    
    return(drop(result))
}

cv.GAMBoost <- function(x=NULL,y,x.linear=NULL,maxstepno=500,family=binomial(),weights=rep(1,length(y)),
                        calc.hat=TRUE,calc.se=TRUE,trace=FALSE,
                        K=10,type="loglik",just.criterion=FALSE,...) 
{
    #   consistency checks
    
    if (is.null(x) && is.null(x.linear)) {
        cat("Neither non-parametric components (x) nor parametric components (x.linear) given. Abort.\n")
        return(NULL)
    }

    all.folds <- split(sample(seq(length(y))), rep(1:K,length=length(y)))

    fit.quality <- matrix(NA,K,maxstepno)

    for (i in 1:K) {
        if (trace) cat("CV-fold",i,":\n")
        omit <- all.folds[[i]]

        fit <- GAMBoost(x=x[-omit,],y[-omit],x.linear=x.linear[-omit,],stepno=maxstepno,family=family,weights=weights[-omit],
                        calc.hat=FALSE,calc.se=FALSE,trace=trace,...)

        prediction <- predict(fit,newdata=x[omit,,drop=FALSE],newdata.linear=x.linear[omit,,drop=FALSE],
                              at.step=1:maxstepno,type="response")

        if (type=="loglik") {
            fit.quality[i,] <- apply(family$dev.resids(matrix(rep(y[omit],maxstepno),length(omit),maxstepno),
                                                        prediction,matrix(rep(weights[omit],maxstepno),length(omit),maxstepno)),2,mean)
        } else {
            if (family$family == "binomial") prediction <- ifelse(prediction > 0.5,1,0)
            fit.quality[i,] <- apply((matrix(rep(y[omit],maxstepno),length(omit),maxstepno) - prediction)^2*weights,2,mean)
        }
    }

    fit.quality.se <- sqrt(apply(fit.quality, 2, var)/K)
    fit.quality <- apply(fit.quality,2,mean)
    #best.fit <- which(fit.quality <= (fit.quality+2*fit.quality.se)[which.min(fit.quality)])[1]
    best.fit <- which.min(fit.quality)
    if (trace) cat("best fit at",best.fit,"\n")

    if (just.criterion) return(list(criterion=fit.quality,se=fit.quality.se,selected=best.fit))
    
    if (trace) cat("final fit:\n")
    return(GAMBoost(x=x,y=y,x.linear=x.linear,stepno=best.fit,family=family,weights=weights,
                    calc.hat=calc.hat,calc.se=calc.se,trace=trace,...))    
}

optimGAMBoostPenalty <- function(x=NULL,y,x.linear=NULL,minstepno=50,maxstepno=200,start.penalty=500,method="AICmin",
                             penalty=100,penalty.linear=100,
                             just.penalty=FALSE,iter.max=10,upper.margin=0.05,trace=TRUE,
                             calc.hat=TRUE,calc.se=TRUE,which.penalty=ifelse(!is.null(x),"smoothness","linear"),...)
{
    #   consistency checks
    
    if (is.null(x) && is.null(x.linear)) {
        cat("Neither non-parametric components (x) nor parametric components (x.linear) given. Abort.\n")
        return(NULL)
    }

    actual.calc.hat <- ifelse(method=="CVmin",FALSE,TRUE)
    actual.calc.se <- ifelse(just.penalty,FALSE,ifelse(method=="CVmin",FALSE,calc.se))
    actual.penalty <- start.penalty
    
    #   default: start from a large penalty and go down, when gone to far use small steps up
    step.up <- 1.2
    step.down <- 0.5
    
    for (i in 1:iter.max) {
        if (trace) cat("iteration",i,": evaluating penalty",actual.penalty,"\n")
        
        if (which.penalty == "smoothness") {
            smoothness.penalty <- actual.penalty
            linear.penalty <- penalty.linear
        } else {
            smoothness.penalty <- penalty
            linear.penalty <- actual.penalty
        }
        
        if (method == "AICmin") {
            actual.res <- GAMBoost(x=x,y,x.linear=x.linear,stepno=maxstepno,penalty=smoothness.penalty,penalty.linear=linear.penalty,
                                   calc.hat=actual.calc.hat,
                                   calc.se=actual.calc.se,trace=trace,...)
            actual.min <- which.min(actual.res$AIC)
        } 
        if (method == "CVmin") {
            actual.res <- cv.GAMBoost(x=x,y,x.linear=x.linear,maxstepno=maxstepno,penalty=smoothness.penalty,penalty.linear=linear.penalty,
                                      just.criterion=TRUE,
                                      calc.hat=actual.calc.hat,calc.se=actual.calc.se,trace=trace,...)
            actual.min <- actual.res$selected
        }
        
        cat("minimum at",actual.min,"\n")
        if (actual.min >= minstepno && actual.min < maxstepno*(1-upper.margin)) break

        #   check whether we are in a scenario where penalty is far to low to start with
        if (i == 1 && actual.min < minstepno) {
            step.up <- 2
            step.down <- 0.8
        }

        if (actual.min < minstepno) {
            actual.penalty <- actual.penalty * step.up
        } else {
            actual.penalty <- actual.penalty * step.down
        }
    }
    
    if (just.penalty) return(actual.penalty)

    if (method == "CVmin") {
        if (which.penalty == "smoothness") {
            smoothness.penalty <- actual.penalty
            linear.penalty <- penalty.linear
        } else {
            smoothness.penalty <- penalty
            linear.penalty <- actual.penalty
        }
                
        actual.res <- GAMBoost(x=x,y,x.linear=x.linear,stepno=maxstepno,penalty=smoothness.penalty,penalty.linear=linear.penalty,
                               calc.hat=calc.hat,calc.se=calc.se,trace=trace,...)
    }
    
    return(actual.res)
}
