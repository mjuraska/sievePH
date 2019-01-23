# 've' returns vaccine efficacy values given parameters for the grid of mark values, v,
# and values for alpha, beta, and gamma
VE <- function(v, alpha, beta, gamma){ 1 - exp(alpha + beta*v + gamma) }

# 'covEstAIPW' returns the estimated covariance matrix of 'phiHat' and 'lambdaHat' in
# the augmented IPW scenario, using Theorem 1 in Juraska and Gilbert (2013, Biometrics)
# 'eventTime' is the observed time, defined as the minimum of failure, censoring, and study times
# 'eventType' is the failure indicator (0 if censored, 1 if failure)
# 'mark' is a data frame (with the same number of rows as the length of 'eventTime') specifying a multivariate mark (a numeric vector for a univariate mark is allowed), with NA for subjects with find=0.
# 'tx' is the treatment group indicator (1 if treatment, 0 if control)
# 'aux' is a numeric vector (a single auxilliary covariate allowed)
# 'auxMiss'
# 'phiHat' is a vector of the alpha and beta estimates
# 'lambdaHat' is the estimate for lambda in the mark density ratio model
# 'gammaHat' is the estimate for gamma obtained in the marginal hazards model
covEstAIPW <- function(eventTime, eventType, mark, tx, aux, auxMiss, phiHat, lambdaHat, gammaHat){
  # convert either a numeric vector or a data frame into a matrix
  mark <- as.matrix(mark)

  n <- length(eventTime)
  isNA <- apply(mark, 1, function(row){ sum(is.na(row)) })
  m.complete <- sum(eventType==1 & isNA==0)
  m.f <- sum(eventType==1)
  eventTime.f <- eventTime[eventType==1]
  eventTime.fM <- matrix(eventTime.f, nrow=n, ncol=m.f, byrow=TRUE)
  eventTime.complete <- eventTime[eventType==1 & isNA==0]
  V.f <- cbind(1,mark[eventType==1,])
  V.complete <- na.omit(V.f)
  na.idx <- attr(V.complete,"na.action")
  tx.complete <- tx[eventType==1 & isNA==0]
  tx.f <- tx[eventType==1]
  aux.complete <- aux[eventType==1 & isNA==0]
  aux.f <- aux[eventType==1]
  auxMiss.f <- auxMiss[eventType==1]
  eventTime.completeM <- matrix(eventTime.complete, nrow=n, ncol=m.complete, byrow=TRUE)
  VV.complete <- apply(V.complete,1,tcrossprod)
  nmark <- NCOL(V.complete)

  g <- function(phi){ exp(drop(V.complete %*% phi)) }
  dG <- function(phi){ t(g(phi) * V.complete) }
  d2G <- function(phi){ array(t(t(VV.complete)*g(phi)),dim=c(nmark,nmark,m.complete)) }
  dGdG <- function(phi){ array(apply(dG(phi),2,tcrossprod),dim=c(nmark,nmark,m.complete)) }

  score1.complete.vect <- function(phi, lambda){
    (-lambda/(1+lambda*(g(phi)-1)) + tx.complete/g(phi)) * t(dG(phi))
  }
  score2.complete.vect <- function(phi, lambda){
    -(g(phi)-1)/(1+lambda*(g(phi)-1))
  }
  aug.mean1 <- function(phi, lambda){
    U <- score1.complete.vect(phi, lambda)
    predicted.vals <- sapply(1:NCOL(U), function(col){
      fit <- lm(U[,col] ~ tx.complete*aux.complete + I(aux.complete^2))
      predict(fit, data.frame(tx.complete=tx.f, aux.complete=aux.f))
    })
    predicted.vals
  }
  aug.mean2 <- function(phi, lambda){
    U <- score2.complete.vect(phi, lambda)
    fit <- lm(U ~ tx.complete*aux.complete + I(aux.complete^2))
    predict(fit, data.frame(tx.complete=tx.f, aux.complete=aux.f))
  }
  score1.vect <- function(phi, lambda){
    vect <- matrix(0, nrow=m.f, ncol=nmark)
    vect[-na.idx,] <- (-lambda/(pi*(1+lambda*(g(phi)-1))) + tx.complete/(pi*g(phi))) * t(dG(phi))
    vect <- vect + (aug.mean1(phi, lambda) * (1-r/pi.all))
    t(vect)
  }
  xi <- function(gamma){ crossprod(eventTime>=eventTime.fM, tx*exp(gamma*tx)) }
  zeta <- function(gamma){ crossprod(eventTime>=eventTime.fM, exp(gamma*tx)) }
  eta <- drop(xi(gammaHat)/zeta(gammaHat))
  score3.vect <- function(gamma){ tx.f-eta }
  l.vect <- function(gamma){
    survprob.vect <- c(1, summary(survfit(Surv(eventTime,eventType)~1))$surv)
    surv.increm <- survprob.vect[-length(survprob.vect)] - survprob.vect[-1]
    eventTime.fMsq <- eventTime.fM[1:m.f,]
    crossprod(eventTime.f>=eventTime.fMsq, surv.increm*(tx.f*exp(gamma*tx.f) - eta*exp(gamma*tx.f))/zeta(gamma))
  }
  score1 <- function(phi, lambda){
    drop(-lambda * dG(phi) %*% (1/(pi*(1+lambda*(g(phi)-1)))) + dG(phi) %*% (tx.complete/(pi*g(phi))) +
           t(aug.mean1(phi, lambda)) %*% (1-r/pi.all))
  }
  score2 <- function(phi, lambda){
    -sum((g(phi)-1)/(pi*(1+lambda*(g(phi)-1)))) + sum(aug.mean2(phi, lambda)*(1-r/pi.all))
  }
  score <- function(phi, lambda){ c(score1(phi,lambda),score2(phi,lambda)) }
  jack11 <- function(phi, lambda){
    d2Gperm <- aperm(d2G(phi), c(3,1,2))
    dGdGperm <- aperm(dGdG(phi), c(3,1,2))
    term1 <- apply(aperm(d2Gperm*(1/(pi*(1+lambda*(g(phi)-1)))), c(2,3,1)),c(1,2),sum)
    term2 <- apply(aperm(dGdGperm*(1/(pi*(1+lambda*(g(phi)-1))^2)), c(2,3,1)),c(1,2),sum)
    term3 <- apply(aperm(d2Gperm*(tx.complete/(pi*g(phi))), c(2,3,1)),c(1,2),sum)
    term4 <- apply(aperm(dGdGperm*(tx.complete/(pi*g(phi)^2)), c(2,3,1)),c(1,2),sum)

    d2U.phi1 <- aperm(d2Gperm*(1/(1+lambda*(g(phi)-1))), c(2,3,1))
    d2U.phi2 <- aperm(dGdGperm*(1/(1+lambda*(g(phi)-1))^2), c(2,3,1))
    d2U.phi3 <- aperm(d2Gperm*(tx.complete/g(phi)), c(2,3,1))
    d2U.phi4 <- aperm(dGdGperm*(tx.complete/g(phi)^2), c(2,3,1))
    d2U.phi <- -lambda*(d2U.phi1 - lambda*d2U.phi2) + d2U.phi3 - d2U.phi4
    predicted.vals.jack <- array(0, dim=c(nmark,nmark,m.f))
    for (i in 1:nmark){
      for (j in 1:nmark){
        resp <- d2U.phi[i,j,]
        fit <- lm(resp ~ tx.complete*aux.complete + I(aux.complete^2))
        predicted.vals.jack[i,j,] <- predict(fit, data.frame(tx.complete=tx.f, aux.complete=aux.f))
      }}
    weighted.predicted.vals.jack <- apply(aperm(aperm(predicted.vals.jack, c(3,1,2))*(1-r/pi.all), c(2,3,1)), c(1,2), sum)
    -lambda*(term1 - lambda*term2) + term3 - term4 + weighted.predicted.vals.jack
  }
  jack21 <- function(phi, lambda){
    d2U.phi.lambda <- -t(dG(phi)) * (1/(1+lambda*(g(phi)-1))^2)
    predicted.vals.jack <- matrix(0,nrow=m.f,ncol=nmark)
    for (i in 1:nmark){
      resp <- d2U.phi.lambda[,i]
      fit <- lm(resp ~ tx.complete*aux.complete + I(aux.complete^2))
      predicted.vals.jack[,i] <- predict(fit, data.frame(tx.complete=tx.f, aux.complete=aux.f))
    }
    weighted.predicted.vals.jack <- colSums(predicted.vals.jack*(1-r/pi.all))
    drop(-dG(phi) %*% (1/(pi*(1+lambda*(g(phi)-1))^2))) + weighted.predicted.vals.jack
  }
  jack22 <- function(phi, lambda){
    d2U.lambda <- ((g(phi)-1)/(1+lambda*(g(phi)-1)))^2
    fit <- lm(d2U.lambda ~ tx.complete*aux.complete + I(aux.complete^2))
    predicted.vals.jack <- predict(fit, data.frame(tx.complete=tx.f, aux.complete=aux.f))
    weighted.predicted.vals.jack <- sum(predicted.vals.jack*(1-r/pi.all))
    sum(((g(phi)-1)^2)/(pi*(1+lambda*(g(phi)-1))^2)) + weighted.predicted.vals.jack
  }
  jack <- function(phi, lambda){
    j21 <- jack21(phi,lambda)
    (cbind(rbind(jack11(phi,lambda),j21),c(j21,jack22(phi,lambda))))/n
  }
  jack33 <- sum(eta*(eta-1))/n

  r <- apply(V.f, 1, function(row){ ifelse(sum(is.na(row))>0,0,1) })
  pi.all <- glm(r ~ tx.f*auxMiss.f, family=binomial)$fitted
  if (!is.null(na.idx)){
    pi <- pi.all[-na.idx]
  } else {
    pi <- pi.all
  }
  if (any(pi<0.005)){ stop("Selection probabilities not bounded away from 0.") }

  p <- mean(eventType==1)
  omega <- drop(score1.vect(phiHat,lambdaHat) %*% (score3.vect(gammaHat) + p*l.vect(gammaHat))/n -
                  sum(score3.vect(gammaHat) + p*l.vect(gammaHat))*apply(score1.vect(phiHat,lambdaHat),1,sum)/(n^2))
  drop(solve(jack(phiHat,lambdaHat))[1:nmark,1:nmark] %*% omega)/(n*jack33)
}

# 'densRatioAIPW' applies the mark density ratio model with missing multivariate marks using
# the inferential procedure defined by augmenting the IPW estimating functions by leveraging
# auxiliary data predictive of the mark.
# The function calculates the mark density ratio and returns a list containing:
#     'coef': estimates for alpha, beta, and lambda
#     'var': the corresponding covariance matrix
#     'jack': the first two rows and columns of the limit estimating function in matrix form
#     'conv': a logical value indicating convergence of the estimating functions
# 'mark' is a data frame representing a multivariate mark variable (a numeric vector for a univariate mark is allowed)
# 'tx' is the treatment group indicator (1 if treatment, 0 if control)
# 'aux' is a numeric vector and it MUST be specified
# 'auxMiss' is a numeric vector and it can be left unspecified
densRatioAIPW <- function(mark, tx, aux, auxMiss=NULL){
  if (missing(aux)){ stop("The auxiliary variable 'aux' for predicting the mean profile score is missing.") }

  # convert either a numeric vector or a data frame into a matrix
  mark <- as.matrix(mark)

  V <- cbind(1,mark)
  V.complete <- na.omit(V)
  z <- tx
  na.idx <- attr(V.complete,"na.action")
  if (!is.null(na.idx)){
    z.complete <- z[-na.idx]
    aux.complete <- aux[-na.idx]
  } else {
    z.complete <- z
    aux.complete <- aux
  }
  nmark <- NCOL(V.complete)
  ninf <- NROW(V.complete)
  VV.complete <- apply(V.complete,1,tcrossprod)

  g <- function(theta){ exp(drop(V.complete %*% theta)) }
  dG <- function(theta){ t(g(theta) * V.complete) }
  d2G <- function(theta){ array(t(t(VV.complete)*g(theta)),dim=c(nmark,nmark,ninf)) }
  dGdG <- function(theta){ array(apply(dG(theta),2,tcrossprod),dim=c(nmark,nmark,ninf)) }

  fscore.i1 <- function(theta, lambda){
    -lambda * t(dG(theta)) * (1/(1+lambda*(g(theta)-1))) + t(dG(theta)) * (z.complete/g(theta))
  }
  fscore.i2 <- function(theta, lambda){
    -(g(theta)-1)/(1+lambda*(g(theta)-1))
  }
  fscore.i <- function(theta, lambda){
    cbind(fscore.i1(theta, lambda), fscore.i2(theta, lambda))
  }
  aug.mean1 <- function(theta, lambda){
    U <- fscore.i1(theta, lambda)
    predicted.vals <- sapply(1:NCOL(U), function(col){
      fit <- lm(U[,col] ~ z.complete*aux.complete + I(aux.complete^2))
      predict(fit, data.frame(z.complete=z, aux.complete=aux))
    })
    return( predicted.vals )
  }
  aug.mean2 <- function(theta, lambda){
    U <- fscore.i2(theta, lambda)
    fit <- lm(U ~ z.complete*aux.complete + I(aux.complete^2))
    predict(fit, data.frame(z.complete=z, aux.complete=aux))
  }

  score1 <- function(theta, lambda){
    drop(-lambda * dG(theta) %*% (1/(pi*(1+lambda*(g(theta)-1)))) + dG(theta) %*% (z.complete/(pi*g(theta))) +
           t(aug.mean1(theta, lambda)) %*% (1-r/pi.all))
  }
  score2 <- function(theta, lambda){
    -sum((g(theta)-1)/(pi*(1+lambda*(g(theta)-1)))) + sum(aug.mean2(theta, lambda)*(1-r/pi.all))
  }
  score <- function(theta, lambda){ c(score1(theta,lambda), score2(theta,lambda)) }
  jack11 <- function(theta, lambda){
    d2Gperm <- aperm(d2G(theta), c(3,1,2))
    dGdGperm <- aperm(dGdG(theta), c(3,1,2))
    term1 <- apply(aperm(d2Gperm*(1/(pi*(1+lambda*(g(theta)-1)))), c(2,3,1)),c(1,2),sum)
    term2 <- apply(aperm(dGdGperm*(1/(pi*(1+lambda*(g(theta)-1))^2)), c(2,3,1)),c(1,2),sum)
    term3 <- apply(aperm(d2Gperm*(z.complete/(pi*g(theta))), c(2,3,1)),c(1,2),sum)
    term4 <- apply(aperm(dGdGperm*(z.complete/(pi*g(theta)^2)), c(2,3,1)),c(1,2),sum)

    d2U.theta1 <- aperm(d2Gperm*(1/(1+lambda*(g(theta)-1))), c(2,3,1))
    d2U.theta2 <- aperm(dGdGperm*(1/(1+lambda*(g(theta)-1))^2), c(2,3,1))
    d2U.theta3 <- aperm(d2Gperm*(z.complete/g(theta)), c(2,3,1))
    d2U.theta4 <- aperm(dGdGperm*(z.complete/g(theta)^2), c(2,3,1))
    d2U.theta <- -lambda*(d2U.theta1 - lambda*d2U.theta2) + d2U.theta3 - d2U.theta4
    predicted.vals.jack <- array(0, dim=c(nmark,nmark,length(z)))
    for (i in 1:nmark){
      for (j in 1:nmark){
        resp <- d2U.theta[i,j,]
        fit <- lm(resp ~ z.complete*aux.complete + I(aux.complete^2))
        predicted.vals.jack[i,j,] <- predict(fit, data.frame(z.complete=z, aux.complete=aux))
      }}
    weighted.predicted.vals.jack <- apply(aperm(aperm(predicted.vals.jack, c(3,1,2))*(1-r/pi.all), c(2,3,1)), c(1,2), sum)

    return( -lambda*(term1 - lambda*term2) + term3 - term4 + weighted.predicted.vals.jack )
  }
  jack21 <- function(theta, lambda){
    d2U.theta.lambda <- -t(dG(theta)) * (1/(1+lambda*(g(theta)-1))^2)
    predicted.vals.jack <- matrix(0,nrow=length(z),ncol=nmark)
    for (i in 1:nmark){
      resp <- d2U.theta.lambda[,i]
      fit <- lm(resp ~ z.complete*aux.complete + I(aux.complete^2))
      predicted.vals.jack[,i] <- predict(fit, data.frame(z.complete=z, aux.complete=aux))
    }
    weighted.predicted.vals.jack <- colSums(predicted.vals.jack*(1-r/pi.all))
    return( drop(-dG(theta) %*% (1/(pi*(1+lambda*(g(theta)-1))^2))) + weighted.predicted.vals.jack )
  }
  jack22 <- function(theta, lambda){
    d2U.lambda <- ((g(theta)-1)/(1+lambda*(g(theta)-1)))^2
    fit <- lm(d2U.lambda ~ z.complete*aux.complete + I(aux.complete^2))
    predicted.vals.jack <- predict(fit, data.frame(z.complete=z, aux.complete=aux))
    weighted.predicted.vals.jack <- sum(predicted.vals.jack*(1-r/pi.all))
    return( sum(((g(theta)-1)^2)/(pi*(1+lambda*(g(theta)-1))^2)) + weighted.predicted.vals.jack )
  }
  jack <- function(theta, lambda){
    j21 <- jack21(theta,lambda)
    cbind(rbind(jack11(theta,lambda),j21),c(j21,jack22(theta,lambda)))
  }

  r <- apply(V, 1, function(row){ ifelse(sum(is.na(row))>0,0,1) })
  if (is.null(auxMiss)) {
    pi.all <- glm(r ~ z, family=binomial)$fitted
  } else {
    pi.all <- glm(r ~ z + auxMiss + z*auxMiss, family=binomial)$fitted
  }
  if (!is.null(na.idx)){
    pi <- pi.all[-na.idx]
  } else {
    pi <- pi.all
  }
  if (any(pi==0)){ stop("Selection probabilities not bounded away from 0.") }

  if (is.null(aux)) {
    param.old <- numeric(nmark+1)
    param.new <- c(numeric(nmark),0.5)
  } else {
    param.old <- densRatioIPW(mark, tx, aux)$coef
    param.new <- param.old + c(numeric(nmark),1e-2)
  }
  while (sum((param.new - param.old)^2)>1e-8){
    param.old <- param.new
    jackInv <- try(solve(jack(param.old[-(nmark+1)],param.old[nmark+1])), silent=TRUE)
    if (class(jackInv)!="try-error"){
      param.new <- param.old - drop(jackInv %*% score(param.old[-(nmark+1)],param.old[nmark+1]))
    }
    if (sum(is.nan(param.new))>0) break
  }
  theta.new <- param.new[-(nmark+1)]
  lambda.new <- param.new[nmark+1]

  Resid <- function(theta, lambda){
    U <- matrix(0,nrow=length(z),ncol=nmark+1)
    U[-na.idx,1:nmark] <- -lambda * t(dG(theta)) * (1/(pi*(1+lambda*(g(theta)-1)))) + t(dG(theta)) * (z.complete/(pi*g(theta)))
    U[,1:nmark] <- U[,1:nmark] + aug.mean1(theta, lambda) * (1-r/pi.all)
    U[-na.idx,nmark+1] <- (g(theta)-1)/(pi*(1+lambda*(g(theta)-1)))
    U[,nmark+1] <- U[,nmark+1] + aug.mean2(theta, lambda) * (1-r/pi.all)
    if (is.null(auxMiss)) {
      S <- (r-pi.all) * cbind(1,z)
      resids <- sapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2])$resid })
    } else {
      S <- (r-pi.all) * cbind(1,z,auxMiss,z*auxMiss)
      resids <- sapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2] + S[,3] + S[,4])$resid })
    }

    crossprod(resids)/ninf
  }

  JackInv <- try(solve(jack(theta.new,lambda.new)), silent=TRUE)
  if (class(JackInv)!="try-error"){
    Var <- ninf * JackInv %*% Resid(theta.new,lambda.new) %*% JackInv
    names(param.new) <- rownames(Var) <- colnames(Var) <- c("alpha",paste("beta",1:(nmark-1),sep=""),"lambda")
  } else {
    Var <- NULL
  }

  list(coef=param.new, var=Var, jack=jack11(theta.new,lambda.new), conv=!(class(jackInv)=="try-error" | class(JackInv)=="try-error"))
}


#' Semiparametric Estimation of Coefficients in a Mark-Specific Proportional Hazards Model
#' with a Missing Multivariate Continuous Mark, with an Inferential Procedure based on 
#' Augmentation of the Inverse Probability Weighted Estimating Functions by Leveraging 
#' Auxiliary Data Predictive of the Mark
#'
#' \code{sievePHaipw} conducts estimation and testing of the multivariate mark-specific hazard ratio
#' model in the competing risks failure time analysis framework for the assessment of mark-specific
#' vaccine efficacy. It accounts for missing multivariate marks by adding an augmented term to the
#' inverse probability weighting (IPW) estimating functions and leveraging auxiliary data predictive
#' of the mark. The method was initially proposed by Robins et al. (1994). Correlations between the
#' mark and auxiliary data are used to "impute' expected profile score vectors for subjects with
#' complete and incomplete mark data.
#' The user can specify whether a one-sided or two-sided hypothesis test is to be performed.
#' 
#' \code{sievePHipw} extends the semiparametric estimation method of Juraska and Gilbert (2013) for continuous mark-
#' specific hazard ratios to accomodate missing at random marks that are univariate or multivariate. The inferential 
#' procedure, initially proposed by Robins et al. (1994), accounts for missing marks by adding an augmented term to the
#' inverse probability weighted (IPW) estimating functions and leveraging auxiliary data predictive of the mark. 
#' Correlations between the mark and auxiliary data are used to "impute" expected profile score vectors for subjects with
#' complete and incomplete mark data and are incorporated into the semiparametric method of maximum profile likelihood 
#' estimation in the treatment-to-placebo mark density ratio model (Qin, 1998). The augmented IPW complete-case estimator 
#' improves efficiency and robustness to mis-specification of the missingness model. \code{sievePHipw} also employs the 
#' ordinary method of maximum partial likelihood estimation of the overall log hazard ratio in the Cox model.
#'
#' @param eventTime a numeric vector specifying the observed right-censored event time
#' @param eventInd a numeric vector indicating the event of interest (1 if event, 0 if right-censored)
#' @param mark either a numeric vector specifying a univariate continuous mark or a data frame specifying a multivariate continuous mark.
#' For subjects with \code{eventInd = 0}, the value(s) in \code{mark} should be set to \code{NA}.
#' @param tx a numeric vector indicating the treatment group (1 if treatment, 0 if placebo)
#' @param aux a numeric vector (a single auxilliary covariate allowed)
#' @param auxMiss 
#'
#' @details
#' \code{sievePHipw} considers data from a randomized placebo-controlled treatment efficacy trial with a time-to-event endpoint.
#' The parameter of interest, the mark-specific hazard ratio, is the ratio (treatment/placebo) of the conditional mark-specific hazard functions.
#' It factors as the product of the mark density ratio (treatment/placebo) and the ordinary marginal hazard function ignoring mark data.
#' The mark density ratio is estimated using the method of Qin (1998) and the Robins et al. (1994) augmented inverse probability weighted (IPW)
#' complete-case estimator, while the marginal hazard ratio is estimated using \code{coxph()} in the \code{survival} package.
#' Both estimators are consistent and asymptotically normal. The asymptotic distribution of the augmented IPW complete-case estimator 
#' is detailed in Juraska and Gilbert (2015).
#'
#' @return An object of class \code{sievePHaipw} which can be processed by
#' \code{\link{summary.sievePHaipw}} to obtain or print a summary of the results. An object of class
#' \code{sievePHaipw} is a list containing the following components:
#' \itemize{
#' \item \code{DRcoef}: a numeric vector of estimates of coefficients \eqn{\phi} in the weight function \eqn{g(v, \phi)} in the density ratio model
#' \item \code{DRlambda}: an estimate of the Lagrange multiplier in the profile score functions for \eqn{\phi} (that arises by profiling out the nuisance parameter)
#' \item \code{DRconverged}: a logical value indicating whether the estimation procedure in the density ratio model converged
#' \item \code{logHR}: an estimate of the marginal log hazard ratio from \code{coxph()} in the \code{survival} package
#' \item \code{cov}: the estimated joint covariance matrix of \code{DRcoef} and \code{logHR}
#' \item \code{coxphFit}: an object returned by the call of \code{coxph()}
#' \item \code{nPlaEvents}: the number of events observed in the placebo group
#' \item \code{nPlaEvents}: the number of events observed in the treatment group
#' \item \code{mark}: the input object
#' \item \code{tx}: the input object
#' }
#' 
#' @references Juraska, M. and Gilbert, P. B. (2013), Mark-specific hazard ratio model with multivariate continuous marks: an application to vaccine efficacy. \emph{Biometrics} 69(2):328-337.
#'
#' Juraska, M., and Gilbert, P. B. (2015). Mark-specific hazard ratio model with missing multivariate marks. \emph{Lifetime data analysis} 22(4): 606-25.
#'
#' Qin, J. (1998), Inferences for case-control and semiparametric two-sample density ratio models. \emph{Biometrika} 85, 619-630.
#'
#' Robins, J. M., Rotnitzky, A., and Zhao, L. P. (1994). Estimation of regression coefficients when some regressors are not always observed. \emph{Journal of the American statistical Association} 89(427): 846-866.
#'
#' @examples
#' n <- 500
#' tx <- rep(0:1, each=n/2)
#' tm <- c(rexp(n/2, 0.1), rexp(n/2, 0.1 * exp(-0.4)))
#' cens <- runif(n, 0, 15)
#' eventTime <- pmin(tm, cens, 3)
#' eventInd <- as.numeric(tm <= pmin(cens, 3))
#' mark1 <- ifelse(eventInd==1, c(rbeta(n/2, 2, 5), rbeta(n/2, 2, 2)), NA)
#' mark2 <- ifelse(eventInd==1, c(rbeta(n/2, 1, 3), rbeta(n/2, 5, 1)), NA)
#'
#' # fit a model with a univariate mark
#' fit <- sievePHaipw(eventTime, eventInd, mark1, tx)
#'
#' # fit a model with a bivariate mark
#' fit <- sievePHaipw(eventTime, eventInd, data.frame(mark1, mark2), tx)
#'
#' @seealso \code{\link{summary.sievePHaipw}} and \code{\link{testIndepTimeMark}}
#'
#' @import survival
#'
#' @export
sievePHaipw <- function(eventTime, eventInd, mark, tx, aux, auxMiss) {
  if (is.numeric(mark)){ mark <- data.frame(mark) }

  nPlaEvents <- sum(eventInd * (1-tx))
  nTxEvents <- sum(eventInd * tx)

  dRatio <- densRatioAIPW(mark[eventInd==1, ], tx[eventInd==1])

  # fit the Cox proportional hazards model to estimate the marginal hazard ratio
  phReg <- coxph(Surv(eventTime, eventInd) ~ tx)
  
  # the estimate of the marginal log hazard ratio
  gammaHat <- phReg$coef
  
  # the output list
  out <- list(DRcoef=NA, DRlambda=NA, DRconverged=dRatio$conv, logHR=gammaHat, cov=NA, coxphFit=phReg, nPlaEvents=nPlaEvents, nTxEvents=nTxEvents, mark=mark, tx=tx)
  
  if (dRatio$conv){

    # a vector of estimates of the density ratio coefficients (alpha, beta1, beta2,..., betak) and the Lagrange multiplier
    thetaHat <- dRatio$coef

    # variance and covariance estimates
    # order of columns in 'dRatio$var': alpha, beta1, beta2,...betak, lambda, where k is the dimension of the mark
    lastComp <- length(thetaHat)
    vthetaHat <- dRatio$var[-lastComp, -lastComp]
    vgammaHat <- drop(phReg$var)
    covThG <- covEstAIPW(eventTime, eventInd, mark, tx, aux, auxMiss, thetaHat[-lastComp], thetaHat[lastComp], gammaHat)

    # covariance matrix for alpha, beta1, beta2,..., betak, gamma
    Sigma <- cbind(rbind(vthetaHat, covThG), c(covThG, vgammaHat))
    colnames(Sigma) <- rownames(Sigma) <- c("alpha", paste0("beta", 1:NCOL(mark)), "gamma")

    out$DRcoef <- thetaHat[-lastComp]
    out$DRlambda <- thetaHat[lastComp]
    out$cov <- Sigma
  }

  class(out) <- "sievePHaipw"
  return(out)

}
