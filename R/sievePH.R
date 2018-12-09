# 've' returns vaccine efficacy values given parameters for the grid of mark values, v,
# and values for alpha, beta, and gamma
VE <- function(v, a, b, g){ 1-exp(a+b*v+g) }

# 'covEst' is an estimator for the covariance matrix of \hat{\phi} and \hat{\lambda} using Theorem 1 in Juraska and Gilbert (2013, Biometrics)
# 'time' is the observed time, defined as the minimum of failure, censoring, and study times
# 'find' is the failure indicator (0 if censored, 1 if failure)
# 'mark' is the mark variable
# 'txind' is the treatment group index (1 if treatment, 0 if control)
# 'phi.hat' is a vector of the alpha and beta estimates
# 'lambda.hat' is the estimate for lambda in the mark density ratio model
# 'gamma.hat' is the estimate for gamma obtained in the marginal hazards model
covEst <- function(time, find, mark, txind, phi.hat, lambda.hat, gamma.hat){
  n <- length(time)
  m <- sum(find)
  time.f <- time[find==1]
  V.f <- cbind(1,mark[find==1])
  txind.f <- txind[find==1]
  time.fM <- matrix(time.f, nrow=n, ncol=m, byrow=TRUE)
  VV <- apply(V.f,1,tcrossprod)
  nmark <- NCOL(V.f)
  
  g <- function(phi){ exp(drop(V.f %*% phi)) }
  dG <- function(phi){ t(g(phi) * V.f) }
  d2G <- function(phi){ array(t(t(VV)*g(phi)),dim=c(nmark,nmark,m)) }
  dGdG <- function(phi){ array(apply(dG(phi),2,tcrossprod),dim=c(nmark,nmark,m)) }
  
  score1.vect <- function(phi, lambda){
    t((-lambda/(1+lambda*(g(phi)-1)) + txind.f/g(phi)) * t(dG(phi)))
  }
  xi <- function(gamma){ crossprod(time>=time.fM, txind*exp(gamma*txind)) }
  zeta <- function(gamma){ crossprod(time>=time.fM, exp(gamma*txind)) }
  eta <- drop(xi(gamma.hat)/zeta(gamma.hat))             
  score3.vect <- function(gamma){ txind.f-eta }
  l.vect <- function(gamma){
    survprob.vect <- c(1, summary(survfit(Surv(time,find)~1), times=sort(time.f))$surv)
    surv.increm <- survprob.vect[-length(survprob.vect)] - survprob.vect[-1]
    time.fMsq <- time.fM[1:m,]
    crossprod(time.f>=time.fMsq, surv.increm*(txind.f*exp(gamma*txind.f) - eta*exp(gamma*txind.f))/zeta(gamma))
  }
  score1 <- function(phi, lambda){
    drop(-lambda * dG(phi) %*% (1/(1+lambda*(g(phi)-1))) + dG(phi) %*% (txind/g(phi)))
  }
  score2 <- function(phi, lambda){
    -sum((g(phi)-1)/(1+lambda*(g(phi)-1)))
  }
  score <- function(phi, lambda){ c(score1(phi,lambda),score2(phi,lambda)) }
  jack11 <- function(phi, lambda){
    d2Gperm <- aperm(d2G(phi), c(3,1,2))
    dGdGperm <- aperm(dGdG(phi), c(3,1,2))
    term1 <- apply(aperm(d2Gperm*(1/(1+lambda*(g(phi)-1))), c(2,3,1)),c(1,2),sum)
    term2 <- apply(aperm(dGdGperm*(1/(1+lambda*(g(phi)-1))^2), c(2,3,1)),c(1,2),sum)
    term3 <- apply(aperm(d2Gperm*(txind.f/g(phi)), c(2,3,1)),c(1,2),sum)
    term4 <- apply(aperm(dGdGperm*(txind.f/g(phi)^2), c(2,3,1)),c(1,2),sum)
    -lambda*(term1 - lambda*term2) + term3 - term4
  }
  jack21 <- function(phi, lambda){
    drop(-dG(phi) %*% (1/(1+lambda*(g(phi)-1))^2))
  }
  jack22 <- function(phi, lambda){
    sum(((g(phi)-1)/(1+lambda*(g(phi)-1)))^2)
  }
  jack <- function(phi, lambda){
    j21 <- jack21(phi,lambda)
    (cbind(rbind(jack11(phi,lambda),j21),c(j21,jack22(phi,lambda))))/n
  }
  jack33 <- sum(eta*(eta-1))/n
  
  p <- mean(find)
  # a vector with 2 components
  omega <- drop(score1.vect(phi.hat,lambda.hat) %*% (score3.vect(gamma.hat) + p*l.vect(gamma.hat))/n - 
                  sum(score3.vect(gamma.hat) + p*l.vect(gamma.hat))*apply(score1.vect(phi.hat,lambda.hat),1,sum)/(n^2))   
  drop(solve(jack(phi.hat,lambda.hat))[1:2,1:2] %*% omega)/(n*jack33)
}

# 'densRatio' calculates the mark density ratio
# 'mark' is the mark variable, which is only observed for cases
# 'trt.id' is the treatment group index (1 if treatment, 0 if control)
densRatio <- function(mark, trt.id){
  V <- cbind(1,mark)
  z <- trt.id
  nmark <- NCOL(V)
  ninf <- NROW(V)
  VV <- apply(V,1,tcrossprod)
  
  g <- function(theta){ exp(drop(V %*% theta)) }
  dG <- function(theta){ t(g(theta) * V) }
  d2G <- function(theta){ array(t(t(VV)*g(theta)),dim=c(nmark,nmark,ninf)) }
  dGdG <- function(theta){ array(apply(dG(theta),2,tcrossprod),dim=c(nmark,nmark,ninf)) }
  
  # profile score functions for the paramter of interest, theta, 
  # and the Lagrange multiplier, lambda
  score1 <- function(theta, lambda){
    drop(-lambda * dG(theta) %*% (1/(1+lambda*(g(theta)-1))) + dG(theta) %*% (z/g(theta)))
  }
  score2 <- function(theta, lambda){
    -sum((g(theta)-1)/(1+lambda*(g(theta)-1)))
  }
  score <- function(theta, lambda){ c(score1(theta,lambda),score2(theta,lambda)) }
  jack11 <- function(theta, lambda){
    d2Gperm <- aperm(d2G(theta), c(3,1,2))
    dGdGperm <- aperm(dGdG(theta), c(3,1,2))
    term1 <- apply(aperm(d2Gperm*(1/(1+lambda*(g(theta)-1))), c(2,3,1)),c(1,2),sum)
    term2 <- apply(aperm(dGdGperm*(1/(1+lambda*(g(theta)-1))^2), c(2,3,1)),c(1,2),sum)
    term3 <- apply(aperm(d2Gperm*(z/g(theta)), c(2,3,1)),c(1,2),sum)
    term4 <- apply(aperm(dGdGperm*(z/g(theta)^2), c(2,3,1)),c(1,2),sum)
    -lambda*(term1 - lambda*term2) + term3 - term4
  }
  jack21 <- function(theta, lambda){
    drop(-dG(theta) %*% (1/(1+lambda*(g(theta)-1))^2))
  }
  jack22 <- function(theta, lambda){
    sum(((g(theta)-1)/(1+lambda*(g(theta)-1)))^2)
  }
  jack <- function(theta, lambda){
    j21 <- jack21(theta,lambda)
    cbind(rbind(jack11(theta,lambda),j21),c(j21,jack22(theta,lambda)))
  }
  
  param.old <- numeric(nmark+1)
  param.new <- c(numeric(nmark),0.5)
  while (sum((param.new - param.old)^2)>1e-8){
    param.old <- param.new
    jackInv <- try(solve(jack(param.old[-(nmark+1)],param.old[nmark+1])), silent=TRUE)
    if (class(jackInv)!="try-error"){
      param.new <- param.old - drop(jackInv %*% score(param.old[-(nmark+1)],
                                                      param.old[nmark+1]))
    }
    if (sum(is.nan(param.new))>0) break
  }
  theta.new <- param.new[-(nmark+1)]
  lambda.new <- param.new[nmark+1]
  
  SigmaHat <- function(theta, lambda){
    L <- -lambda * t(dG(theta)) * (1/(1+lambda*(g(theta)-1))) + t(dG(theta)) * (z/g(theta))
    L <- cbind(L, (g(theta)-1)/(1+lambda*(g(theta)-1)))
    crossprod(L)/ninf
  }
  
  JackInv <- try(solve(jack(theta.new,lambda.new)), silent=TRUE)
  if (class(JackInv)!="try-error"){
    Var <- ninf * JackInv %*% SigmaHat(theta.new,lambda.new) %*% JackInv
    names(param.new) <- rownames(Var) <- colnames(Var) <- c("alpha",
                                                            paste("beta",1:(nmark-1),sep=""),"lambda")
  } else {
    Var <- NULL
  }
  
  list(coef=param.new, var=Var, jack=jack11(theta.new,lambda.new),
       conv=!(class(jackInv)=="try-error" | class(JackInv)=="try-error"))
}  


#' Semiparametric Efficient Estimation and Testing for a Mark-Specific Proportional Hazards 
#' Model with Multivariate Continuous Marks
#' 
#' \code{sievePH} conducts estimation and testing of the multivariate mark-specific hazard 
#' ratio model in the competing risks failure time analysis framework for the assessment of 
#' mark-specific vaccine efficacy. It improves efficiency by employing the semiparametric 
#' method of maximum profile likelihood estimation in the vaccine-to-placebo mark density 
#' ratio model. The method employs a complete-cases analysis of the mark in failures. 
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
sievePH <- function(time, failInd, mark, txInd, oneSided) {
  X <- time
  d <- failInd
  V <- mark
  Z <- txInd
  
  dRatio <- densRatio(V[d==1],Z[d==1])
  
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
    covThG <- covEst(X,d,V,Z,thetaHat[1:2],thetaHat[3],gammaHat)
    
    # vaccine efficacy estimate
    ve <- VE(V,thetaHat[1],thetaHat[2],gammaHat)
    
    # covariance matrix for alpha, beta1, covThg
    Sigma <- cbind(rbind(vthetaHat,covThG), c(covThG,vgammaHat)) 
    # colnames(Sigma) <- rownames(Sigma) <- c("alpha", "beta1", "gamma")
    
    return( list(mark = V, txInd = Z, nInf0 = nInf0, nInf1 = nInf1, alphaHat=thetaHat[1], betaHat=thetaHat[2], lambdaHat = thetaHat[3], gammaHat = gammaHat, 
                 ve, Sigma))
  }
  
}
  