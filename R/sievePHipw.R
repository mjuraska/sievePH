# 've' returns vaccine efficacy values given parameters for the grid of mark values, v,
# and values for alpha, beta, and gamma
VE <- function(v, alpha, beta, gamma){ 1 - exp(alpha + beta*v + gamma) }

# 'covEstIPW' returns the estimated covariance matrix of 'phiHat' and 'lambdaHat'
# in the IPW scenario, using Theorem 1 in Juraska and Gilbert (2013, Biometrics)
# 'eventTime' is the observed time, defined as the minimum of failure, censoring, and study times
# 'eventType' is the failure indicator (0 if censored, 1 if failure)
# 'mark' is a data frame (with the same number of rows as the length of 'eventTime') specifying a multivariate mark (a numeric vector for a univariate mark is allowed), with NA for subjects with find=0.
# 'tx' is the treatment group indicator (1 if treatment, 0 if control)
# 'auxMiss'
# 'phiHat' is a vector of the alpha and beta estimates
# 'lambdaHat' is the estimate for lambda in the mark density ratio model
# 'gammaHat' is the estimate for gamma obtained in the marginal hazards model
covEstIPW <- function(eventTime, eventType, mark, tx, auxMiss, phiHat, lambdaHat, gammaHat){
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
  auxMiss.f <- auxMiss[eventType==1]
  eventTime.completeM <- matrix(eventTime.complete, nrow=n, ncol=m.complete, byrow=TRUE)
  VV.complete <- apply(V.complete,1,tcrossprod)
  nmark <- NCOL(V.complete)

  g <- function(phi){ exp(drop(V.complete %*% phi)) }
  dG <- function(phi){ t(g(phi) * V.complete) }
  d2G <- function(phi){ array(t(t(VV.complete)*g(phi)),dim=c(nmark,nmark,m.complete)) }
  dGdG <- function(phi){ array(apply(dG(phi),2,tcrossprod),dim=c(nmark,nmark,m.complete)) }

  score1.vect <- function(phi, lambda){
    vect <- matrix(0, nrow=nmark, ncol=m.f)
    vect[,-na.idx] <- t((-lambda/(pi*(1+lambda*(g(phi)-1))) + tx.complete/(pi*g(phi))) * t(dG(phi)))
    vect
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
    drop(-lambda * dG(phi) %*% (1/(pi*(1+lambda*(g(phi)-1)))) + dG(phi) %*% (tx.complete/(pi*g(phi))))
  }
  score2 <- function(phi, lambda){
    -sum((g(phi)-1)/(pi*(1+lambda*(g(phi)-1))))
  }
  score <- function(phi, lambda){ c(score1(phi,lambda),score2(phi,lambda)) }
  jack11 <- function(phi, lambda){
    d2Gperm <- aperm(d2G(phi), c(3,1,2))
    dGdGperm <- aperm(dGdG(phi), c(3,1,2))
    term1 <- apply(aperm(d2Gperm*(1/(pi*(1+lambda*(g(phi)-1)))), c(2,3,1)),c(1,2),sum)
    term2 <- apply(aperm(dGdGperm*(1/(pi*(1+lambda*(g(phi)-1))^2)), c(2,3,1)),c(1,2),sum)
    term3 <- apply(aperm(d2Gperm*(tx.complete/(pi*g(phi))), c(2,3,1)),c(1,2),sum)
    term4 <- apply(aperm(dGdGperm*(tx.complete/(pi*g(phi)^2)), c(2,3,1)),c(1,2),sum)
    -lambda*(term1 - lambda*term2) + term3 - term4
  }
  jack21 <- function(phi, lambda){
    drop(-dG(phi) %*% (1/(pi*(1+lambda*(g(phi)-1))^2)))
  }
  jack22 <- function(phi, lambda){
    sum(((g(phi)-1)^2)/(pi*(1+lambda*(g(phi)-1))^2))
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

# 'densRatioIPW' applies the mark density ratio model with missing multivariate marks using
# the inferential procedure defined by inverse probability weighting (IPW) of complete cases.
# The function calculates the mark density ratio and returns a list containing:
#     'coef': estimates for alpha, beta, and lambda
#     'var': the corresponding covariance matrix
#     'jack': the first two rows and columns of the limit estimating function in matrix form
#     'conv': a logical value indicating convergence of the estimating functions
# 'mark' is a data frame representing a multivariate mark variable (a numeric vector for a univariate mark is allowed)
# 'tx' is the treatment group indicator (1 if treatment, 0 if control)
# 'aux' is a numeric vector (a single auxilliary covariate allowed)
densRatioIPW <- function(mark, tx, aux=NULL){
  # convert either a numeric vector or a data frame into a matrix
  mark <- as.matrix(mark)

  V <- cbind(1,mark)
  V.complete <- na.omit(V)
  z <- tx
  na.idx <- attr(V.complete,"na.action")
  if (!is.null(na.idx)){
    z.complete <- z[-na.idx]
  } else {
    z.complete <- z
  }
  nmark <- NCOL(V.complete)
  ninf <- NROW(V.complete)
  VV.complete <- apply(V.complete,1,tcrossprod)

  g <- function(theta){ exp(drop(V.complete %*% theta)) }
  dG <- function(theta){ t(g(theta) * V.complete) }
  d2G <- function(theta){ array(t(t(VV.complete)*g(theta)),dim=c(nmark,nmark,ninf)) }
  dGdG <- function(theta){ array(apply(dG(theta),2,tcrossprod),dim=c(nmark,nmark,ninf)) }

  score1 <- function(theta, lambda){
    drop(-lambda * dG(theta) %*% (1/(pi*(1+lambda*(g(theta)-1)))) + dG(theta) %*% (z.complete/(pi*g(theta))))
  }
  score2 <- function(theta, lambda){
    -sum((g(theta)-1)/(pi*(1+lambda*(g(theta)-1))))
  }
  score <- function(theta, lambda){ c(score1(theta,lambda),score2(theta,lambda)) }
  jack11 <- function(theta, lambda){
    d2Gperm <- aperm(d2G(theta), c(3,1,2))
    dGdGperm <- aperm(dGdG(theta), c(3,1,2))
    term1 <- apply(aperm(d2Gperm*(1/(pi*(1+lambda*(g(theta)-1)))), c(2,3,1)),c(1,2),sum)
    term2 <- apply(aperm(dGdGperm*(1/(pi*(1+lambda*(g(theta)-1))^2)), c(2,3,1)),c(1,2),sum)
    term3 <- apply(aperm(d2Gperm*(z.complete/(pi*g(theta))), c(2,3,1)),c(1,2),sum)
    term4 <- apply(aperm(dGdGperm*(z.complete/(pi*g(theta)^2)), c(2,3,1)),c(1,2),sum)
    -lambda*(term1 - lambda*term2) + term3 - term4
  }
  jack21 <- function(theta, lambda){
    drop(-dG(theta) %*% (1/(pi*(1+lambda*(g(theta)-1))^2)))
  }
  jack22 <- function(theta, lambda){
    sum(((g(theta)-1)^2)/(pi*(1+lambda*(g(theta)-1))^2))
  }
  jack <- function(theta, lambda){
    j21 <- jack21(theta,lambda)
    cbind(rbind(jack11(theta,lambda),j21),c(j21,jack22(theta,lambda)))
  }

  r <- apply(V, 1, function(row){ ifelse(sum(is.na(row))>0,0,1) })
  if (is.null(aux)){
    pi.all <- glm(r ~ z, family=binomial)$fitted
  } else {
    pi.all <- glm(r ~ z + aux + z*aux, family=binomial)$fitted
  }

  if (!is.null(na.idx)){
    pi <- pi.all[-na.idx]
  } else {
    pi <- pi.all
  }
  if (any(pi<0.005)){ stop("Selection probabilities not bounded away from 0.") }

  param.old <- numeric(nmark+1)
  param.new <- c(numeric(nmark),0.5)
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
    U[-na.idx,nmark+1] <- (g(theta)-1)/(pi*(1+lambda*(g(theta)-1)))

    if (is.null(aux)) {
      S <- (r-pi.all) * cbind(1,z)
      resids <- lapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2])$resid })
    } else {
      S <- (r-pi.all) * cbind(1,z,aux,z*aux)
      resids <- lapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2] + S[,3] + S[,4])$resid })
    }

    Resids <- do.call("cbind",resids)
    crossprod(Resids)/ninf
  }

  JackInv <- try(solve(jack(theta.new,lambda.new)), silent=TRUE)
  if (class(JackInv)!="try-error"){
    Var <- ninf * JackInv %*% Resid(theta.new,lambda.new) %*% JackInv
    names(param.new) <- rownames(Var) <- colnames(Var) <- c("alpha",paste("beta",1:(nmark-1),sep=""),"lambda")
  } else {
    Var <- NULL
  }

  list(coef=param.new, var=Var, jack=jack11(theta.new,lambda.new), probs=pi, conv=!(class(jackInv)=="try-error" | class(JackInv)=="try-error"))
}


#' Semiparametric Estimation of Coefficients in a Mark-Specific Proportional Hazards Model
#' with a Multivariate Continuous Mark, with Missing Values Inferred using Inverse 
#' Probability Weighting of Complete Cases
#'
#' \code{sievePHipw} extends the semiparametric estimation method of Juraska and Gilbert (2013) for continuous mark-
#' specific hazard ratios to accomodate missing at random marks that are univariate or multivariate. The inferential 
#' procedure is based on inverse probability weighting (IPW) of the complete cases, originally proposed by 
#' Horvitz and Thompson (1952), where the complete cases are weighted by the inverse of their estimated probabilities 
#' given their event status, treatment status, and auxiliary covariates. This information is incorporated into weighted 
#' profile score functions in the semiparametric method of maximum profile likelihood estimation in the treatment-to-placebo 
#' mark density ratio model (Qin, 1998). \code{sievePHipw} also employs the ordinary method of maximum partial likelihood 
#' estimation of the overall log hazard ratio in the Cox model.
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
#' The mark density ratio is estimated using the method of Qin (1998) and the Horvitz and Thompson (1952) inverse probability weighted (IPW)
#' complete-case estimator, while the marginal hazard ratio is estimated using \code{coxph()} in the \code{survival} package.
#' Both estimators are consistent and asymptotically normal. The asymptotic distribution of the IPW complete-case estimator is detailed in Juraska and Gilbert (2015).
#' 
#' @return An object of class \code{sievePHipw} which can be processed by
#' \code{\link{summary.sievePHipw}} to obtain or print a summary of the results. An object of class
#' \code{sievePHipw} is a list containing the following components:
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
#' @references Horvitz, D. G., and Thompson, D. J. (1952). A generalization of sampling without replacement from a finite universe. \emph{Journal of the American statistical Association} 47(260): 663-685.
#' 
#' Juraska, M. and Gilbert, P. B. (2013), Mark-specific hazard ratio model with multivariate continuous marks: an application to vaccine efficacy. \emph{Biometrics} 69(2):328-337.
#' 
#' Juraska, M., and Gilbert, P. B. (2015). Mark-specific hazard ratio model with missing multivariate marks. \emph{Lifetime data analysis} 22(4): 606-25.
#' 
#' Qin, J. (1998), Inferences for case-control and semiparametric two-sample density ratio models. \emph{Biometrika} 85, 619-630.
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
#' fit <- sievePHipw(eventTime, eventInd, mark1, tx)
#'
#' # fit a model with a bivariate mark
#' fit <- sievePHipw(eventTime, eventInd, data.frame(mark1, mark2), tx)
#'
#' @seealso \code{\link{summary.sievePHipw}} and \code{\link{testIndepTimeMark}}
#'
#' @import survival
#'
#' @export
sievePHipw <- function(eventTime, eventInd, mark, tx, aux = NULL, auxMiss = NULL) {
  if (is.numeric(mark)){ mark <- data.frame(mark) }
  
  nPlaEvents <- sum(eventInd * (1-tx))
  nTxEvents <- sum(eventInd * tx)

  dRatio <- densRatioIPW(mark[eventInd==1, ], tx[eventInd==1], aux)
  
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
    covThG <- covEstIPW(eventTime, eventInd, mark, tx, auxMiss, thetaHat[-lastComp], thetaHat[lastComp], gammaHat)

    # covariance matrix for alpha, beta1, beta2,..., betak, gamma
    Sigma <- cbind(rbind(vthetaHat, covThG), c(covThG, vgammaHat))
    colnames(Sigma) <- rownames(Sigma) <- c("alpha", paste0("beta", 1:NCOL(mark)), "gamma")

    out$DRcoef <- thetaHat[-lastComp]
    out$DRlambda <- thetaHat[lastComp]
    out$cov <- Sigma
  }

  class(out) <- "sievePHipw"
  return(out)
}
