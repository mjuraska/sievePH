# 've' returns vaccine efficacy values given parameters for the grid of mark values, v,
# and values for alpha, beta, and gamma
VE <- function(v, alpha, beta, gamma){ 1 - exp(alpha + beta*v + gamma) }

# 'covEstIPW' returns the estimated covariance matrix of 'phiHat' and 'lambdaHat' 
# in the IPW scenario, using Theorem 1 in Juraska and Gilbert (2013, Biometrics)
# 'eventTime' is the observed time, defined as the minimum of failure, censoring, and study times
# 'eventType' is the failure indicator (0 if censored, 1 if failure)
# 'mark' is the mark variable
# 'tx' is the treatment group indicator (1 if treatment, 0 if control)
# 'auxMiss'
# 'phiHat' is a vector of the alpha and beta estimates
# 'lambdaHat' is the estimate for lambda in the mark density ratio model
# 'gammaHat' is the estimate for gamma obtained in the marginal hazards model
covEstIPW <- function(eventTime, eventType, mark, tx, auxMiss, phiHat, lambdaHat, gammaHat){
  n <- length(eventTime)
  m.complete <- sum(eventType==1 & !is.na(mark))
  m.f <- sum(eventType==1)
  eventTime.f <- eventTime[eventType==1]
  eventTime.fM <- matrix(eventTime.f, nrow=n, ncol=m.f, byrow=TRUE)
  eventTime.complete <- eventTime[eventType==1 & !is.na(mark)]
  V.f <- cbind(1,mark[eventType==1])
  V.complete <- na.omit(V.f)
  na.idx <- attr(V.complete,"na.action")	
  tx.complete <- tx[eventType==1 & !is.na(mark)]
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
                  sum(score3.vect(gammaHat) + p*l.vect(gammaHat))*apply(score1.vect(phiHat,lambdaHat),1,sum)/(n^2))   # a vector with 2 components
  drop(solve(jack(phiHat,lambdaHat))[1:2,1:2] %*% omega)/(n*jack33)
}

# 'densRatioIPW' applies the mark density ratio model with missing multivariate marks using 
# the inferential procedure defined by inverse probability weighting (IPW) of complete cases.
# The function calculates the mark density ratio and returns a list containing:
#     'coef': estimates for alpha, beta, and lambda 
#     'var': the corresponding covariance matrix 
#     'jack': the first two rows and columns of the limit estimating function in matrix form
#     'conv': a logical value indicating convergence of the estimating functions
# 'mark' is the mark variable, which is only observed for cases
# 'tx' is the treatment group indicator (1 if treatment, 0 if control)
# 'aux'
densRatioIPW <- function(mark, tx, aux){
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
  # pi.all <- glm(r ~ z + aux + z*aux, family=binomial)$fitted
  pi.all <- glm(r ~ z, family=binomial)$fitted
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


#' Semiparametric Efficient Estimation and Testing for a Mark-Specific Proportional Hazards 
#' Model with Missing Multivariate Marks using Inverse Probability Weighting of Complete Cases
#' 
#' \code{sievePH.IPW} conducts estimation and testing of the multivariate mark-specific hazard 
#' ratio model in the competing risks failure time analysis framework for the assessment of 
#' mark-specific vaccine efficacy. It accounts for missing multivariate marks using the 
#' inferential procedure based on inverse probability weighting (IPW) of the complete cases, 
#' originally proposed by Horvitz and Thompson (1952), where the complete cases are weighted 
#' by the inverse of the probabilities or their estimates.
#' The user can specify whether a one-sided or two-sided hypothesis test is to be performed.
#' 
#' @param eventTime a numeric vector specifying the observed time, defined as the minimum of the 
#' event, censoring, and study time.
#' @param eventType a binary vector indicating the event status (1 if failure, 0 if censored)
#' @param mark a numeric vector of the values of the mark variable, observed only in cases
#' @param tx a binary vector indicating the treatment group (1 if treatment, 0 if control)
#' @param aux
#' @param auxMiss
#' 
#' @details 
#' 
#' @return \code{sievePH.IPW} returns an object of class "sievePH.IPW" which can be processed by 
#' \code{\link{summary.sievePH.IPW}} to obtain or print a summary of the results. An object of class
#' "sievePH.IPW" is a list containing the following components:
#' \item{alphaHat}{the estimate for the \eqn{alpha} parameter in the density ratio model} 
#' \item{betaHat}{the estimate for the \eqn{beta} parameter in the density ratio model}
#' \item{gammaHat}{the estimate for the \eqn{gamma} parameter in the density ratio model}
#' \item{lambdaHat}{the estimate for \eqn{lambda}, the Lagrange multiplier utilized in 
#' the estimation procedure}
#' \item{cov}{the covariance matrix for \eqn{alpha}, \eqn{beta}, and \eqn{gamma}}
#' \item{ve}{a numeric vector of estimates of vaccine efficacy}
#' \item{mark}{a numeric vector of the values of the mark variable, observed only in cases}
#' \item{tx}{a binary vector indicating the treatment group (1 if treatment, 0 if control)}
#' \item{nEvents0}{the number of events in the placebo group}
#' \item{nEvents1}{the number of events in the vaccine group}
#' \item{coxModel}{the fitted cox regression model for the marginal hazard ratio}
#' 
#' @examples 
#' 
#' @import survival
#' 
#' @export
sievePH.IPW <- function(eventTime, eventType, mark, tx, aux = NULL, auxMiss = NULL) {

  X <- eventTime
  d <- eventType
  V <- mark
  Z <- tx
  nEvents0 <- sum(d*(1-Z))                             # number of failures in placebo group
  nEvents1 <- sum(d*Z)                                 # number of failures in vaccine group
  
  dRatio <- densRatioIPW(V[d==1],Z[d==1], aux)
  
  if (dRatio$conv){
    
    # Cox proportional hazards model for marginal hazard ratio
    phReg <- coxph(Surv(X,d)~Z)
    
    # parameter estimates
    thetaHat <- dRatio$coef
    gammaHat <- phReg$coef
    
    # variance and covariance estimates
    # order of columns: alpha, beta1, beta2,...betak, lambda, where k is number of marks
    vthetaHat <- dRatio$var[1:2,1:2]
    vgammaHat <- drop(phReg$var) 
    covThG <- covEstIPW(X,d,V,Z,auxMiss,thetaHat[1:2],thetaHat[3],gammaHat)
    
    # vaccine efficacy estimate
    ve <- VE(V,thetaHat[1],thetaHat[2],gammaHat)
    
    # covariance matrix for alpha, beta1, gamma
    Sigma <- cbind(rbind(vthetaHat,covThG), c(covThG,vgammaHat)) 
    colnames(Sigma) <- rownames(Sigma) <- c("alpha", "beta1", "gamma")
    
    result <- list(mark = V, tx = Z, nEvents0 = nEvents0, nEvents1 = nEvents1, alphaHat=thetaHat[1], betaHat=thetaHat[2], lambdaHat = thetaHat[3], gammaHat = gammaHat, 
                   ve = ve, cov = Sigma, coxModel = phReg)
  } else {
    result <- list(mark = V, tx = Z, nEvents0 = nEvents0, nEvents1 = nEvents1)
  }
  
  class(result) <- "sievePH.IPW"
  return(result)
}
