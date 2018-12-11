# 've' returns vaccine efficacy values given parameters for the grid of mark values, v,
# and values for alpha, beta, and gamma
VE <- function(v, a, b, g){ 1-exp(a+b*v+g) }

# 'covEstIPW' is an estimator for \cov(\hat{\phi}_{ipw},\hat{\lambda}) based on Theorem 1 of JG (2013)
# 'time' is the observed time, defined as the minimum of failure, censoring, and study times
# 'find' is the failure indicator (0 if censored, 1 if failure)
# 'mark' is the mark variable
# 'txind' is the treatment group index (1 if treatment, 0 if control)
# 'aux.miss' 
# 'phi.hat' is a vector of the alpha and beta estimates
# 'lambda.hat' is the estimate for lambda in the mark density ratio model
# 'gamma.hat' is the estimate for gamma obtained in the marginal hazards model
covEstIPW <- function(time, find, mark, txind, aux.miss, phi.hat, lambda.hat, gamma.hat){
  n <- length(time)
  m.complete <- sum(find==1 & !is.na(mark))
  m.f <- sum(find==1)
  time.f <- time[find==1]
  time.fM <- matrix(time.f, nrow=n, ncol=m.f, byrow=TRUE)
  time.complete <- time[find==1 & !is.na(mark)]
  V.f <- cbind(1,mark[find==1])
  V.complete <- na.omit(V.f)
  na.idx <- attr(V.complete,"na.action")	
  txind.complete <- txind[find==1 & !is.na(mark)]
  txind.f <- txind[find==1]
  aux.miss.f <- aux.miss[find==1]
  time.completeM <- matrix(time.complete, nrow=n, ncol=m.complete, byrow=TRUE)
  VV.complete <- apply(V.complete,1,tcrossprod)
  nmark <- NCOL(V.complete)
  
  g <- function(phi){ exp(drop(V.complete %*% phi)) }
  dG <- function(phi){ t(g(phi) * V.complete) }
  d2G <- function(phi){ array(t(t(VV.complete)*g(phi)),dim=c(nmark,nmark,m.complete)) }
  dGdG <- function(phi){ array(apply(dG(phi),2,tcrossprod),dim=c(nmark,nmark,m.complete)) }
  
  score1.vect <- function(phi, lambda){
    vect <- matrix(0, nrow=nmark, ncol=m.f)
    vect[,-na.idx] <- t((-lambda/(pi*(1+lambda*(g(phi)-1))) + txind.complete/(pi*g(phi))) * t(dG(phi)))
    vect
  }
  xi <- function(gamma){ crossprod(time>=time.fM, txind*exp(gamma*txind)) }
  zeta <- function(gamma){ crossprod(time>=time.fM, exp(gamma*txind)) }
  eta <- drop(xi(gamma.hat)/zeta(gamma.hat))             
  score3.vect <- function(gamma){ txind.f-eta }
  l.vect <- function(gamma){
    survprob.vect <- c(1, summary(survfit(Surv(time,find)~1))$surv)
    surv.increm <- survprob.vect[-length(survprob.vect)] - survprob.vect[-1]
    time.fMsq <- time.fM[1:m.f,]
    crossprod(time.f>=time.fMsq, surv.increm*(txind.f*exp(gamma*txind.f) - eta*exp(gamma*txind.f))/zeta(gamma))
  }
  score1 <- function(phi, lambda){
    drop(-lambda * dG(phi) %*% (1/(pi*(1+lambda*(g(phi)-1)))) + dG(phi) %*% (txind.complete/(pi*g(phi))))
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
    term3 <- apply(aperm(d2Gperm*(txind.complete/(pi*g(phi))), c(2,3,1)),c(1,2),sum)
    term4 <- apply(aperm(dGdGperm*(txind.complete/(pi*g(phi)^2)), c(2,3,1)),c(1,2),sum)
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
  pi.all <- glm(r ~ txind.f*aux.miss.f, family=binomial)$fitted
  if (!is.null(na.idx)){
    pi <- pi.all[-na.idx]
  } else {
    pi <- pi.all
  }
  if (any(pi<0.005)){ stop("Selection probabilities not bounded away from 0.") }
  
  p <- mean(find==1)
  omega <- drop(score1.vect(phi.hat,lambda.hat) %*% (score3.vect(gamma.hat) + p*l.vect(gamma.hat))/n - 
                  sum(score3.vect(gamma.hat) + p*l.vect(gamma.hat))*apply(score1.vect(phi.hat,lambda.hat),1,sum)/(n^2))   # a vector with 2 components
  drop(solve(jack(phi.hat,lambda.hat))[1:2,1:2] %*% omega)/(n*jack33)
}

# 'densRatioIPW' applies the mark density ratio model with missing multivariate marks using 
# the inferential procedure defined by inverse probability weighting (IPW) of complete cases
# 'mark' is the mark variable, which is only observed for cases
# 'trt.id' is the treatment group index (1 if treatment, 0 if control)
# 'aux' 
densRatioIPW <- function(mark, trt.id, aux){
  V <- cbind(1,mark)
  V.complete <- na.omit(V)
  z <- trt.id
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
    # S <- (r-pi.all) * cbind(1,z,aux,z*aux)
    S <- (r-pi.all) * cbind(1,z)
    # resids <- lapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2] + S[,3] + S[,4])$resid })
    resids <- lapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2])$resid })
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
#' \code{sieveIPW} conducts estimation and testing of the multivariate mark-specific hazard 
#' ratio model in the competing risks failure time analysis framework for the assessment of 
#' mark-specific vaccine efficacy. It accounts for missing multivariate marks using the 
#' inferential procude based on inverse probability weighting (IPW) of the complete cases, 
#' originally proposed by Horvitz and Thompson (1952), where the complete cases are weighted 
#' by the inverse of the probabilities or their estimates.
#' The user can specify whether a one-sided or two-sided hypothesis test is to be performed.
#' 
#' @param time a numeric vector specifying the observed time, defined as the minimum of the failure, censoring, and study time.
#' @param failInd an binary vector indicating the failure status (1 if failure, 0 if censored)
#' @param mark a numeric vector specifying the values of the mark variable, observed only in cases
#' @param txInd a binary vector indicating the treatment group (1 if treatment, 0 if control)
#' @param oneSided a logical value indicating if the the sieve test of H0:HR(v)=HR is to be one-sided (TRUE) or two-sided (FALSE). 
#' 
#' @details 
#' 
#' @return Returns an output list with the following components:
#' \itemize{
#'   \item \code{alphaHat}: 
#' }
#' 
#' @examples 
#' 
#' @import survival
#' 
#' @export
sieveIPW <- function(time, failInd, mark, txInd, oneSided) {

  X <- time
  d <- failInd
  V <- mark
  Z <- txInd
  nFail0 <- sum(d*(1-Z))                             # number of failures in placebo group
  nFail1 <- sum(d*Z)                                 # number of failures in vaccine group
  
  dRatio <- densRatioIPW(V[d==1],Z[d==1])
  
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
    covThG <- covEstIPW(X,d,V,Z,thetaHat[1:2],thetaHat[3],gammaHat)
    
    # vaccine efficacy estimate
    ve <- VE(V,thetaHat[1],thetaHat[2],gammaHat)
    
    # covariance matrix for alpha, beta1, gamma
    Sigma <- cbind(rbind(vthetaHat,covThG), c(covThG,vgammaHat)) 
    colnames(Sigma) <- rownames(Sigma) <- c("alpha", "beta1", "gamma")
    
    return( list(mark = V, txInd = Z, nFail0 = nFail0, nFail1 = nFail1, alphaHat=thetaHat[1], betaHat=thetaHat[2], lambdaHat = thetaHat[3], gammaHat = gammaHat, 
                 ve = ve, cov = Sigma))
  } else {
    return( list(mark = V, txInd = Z, nFail0 = nFail0, nFail1 = nFail1))
  }
  
}
