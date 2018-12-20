# 've' returns vaccine efficacy values given parameters for the grid of mark values, v,
# and values for alpha, beta, and gamma
VE <- function(v, alpha, beta, gamma){ 1 - exp(alpha + beta*v + gamma) }

# 'covEst' returns the estimated covariance matrix of 'phiHat' and 'lambdaHat' using 
# Theorem 1 in Juraska and Gilbert (2013, Biometrics)
# 'time' is the observed time, defined as the minimum of failure, censoring, and study times
# 'find' is the failure indicator (0 if censored, 1 if failure)
# 'mark' is the mark variable
# 'txInd' is the treatment group indicator (1 if treatment, 0 if control)
# 'phiHat' is a vector of the alpha and beta estimates
# 'lambdaHat' is the estimate for lambda in the mark density ratio model
# 'gammaHat' is the estimate for gamma obtained in the marginal hazards model
covEst <- function(time, find, mark, txInd, phiHat, lambdaHat, gammaHat){
  n <- length(time)
  m <- sum(find)
  time.f <- time[find==1]
  V.f <- cbind(1,mark[find==1])
  txInd.f <- txInd[find==1]
  time.fM <- matrix(time.f, nrow=n, ncol=m, byrow=TRUE)
  VV <- apply(V.f,1,tcrossprod)
  nmark <- NCOL(V.f)

  g <- function(phi){ exp(drop(V.f %*% phi)) }
  dG <- function(phi){ t(g(phi) * V.f) }
  d2G <- function(phi){ array(t(t(VV)*g(phi)),dim=c(nmark,nmark,m)) }
  dGdG <- function(phi){ array(apply(dG(phi),2,tcrossprod),dim=c(nmark,nmark,m)) }

  score1.vect <- function(phi, lambda){
    t((-lambda/(1+lambda*(g(phi)-1)) + txInd.f/g(phi)) * t(dG(phi)))
  }
  xi <- function(gamma){ crossprod(time>=time.fM, txInd*exp(gamma*txInd)) }
  zeta <- function(gamma){ crossprod(time>=time.fM, exp(gamma*txInd)) }
  eta <- drop(xi(gammaHat)/zeta(gammaHat))             
  score3.vect <- function(gamma){ txInd.f-eta }
  l.vect <- function(gamma){
    survprob.vect <- c(1, summary(survfit(Surv(time,find)~1), times=sort(time.f))$surv)
    surv.increm <- survprob.vect[-length(survprob.vect)] - survprob.vect[-1]
    time.fMsq <- time.fM[1:m,]
    crossprod(time.f>=time.fMsq, surv.increm*(txInd.f*exp(gamma*txInd.f) - eta*exp(gamma*txInd.f))/zeta(gamma))
  }
  score1 <- function(phi, lambda){
    drop(-lambda * dG(phi) %*% (1/(1+lambda*(g(phi)-1))) + dG(phi) %*% (txInd/g(phi)))
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
    term3 <- apply(aperm(d2Gperm*(txInd.f/g(phi)), c(2,3,1)),c(1,2),sum)
    term4 <- apply(aperm(dGdGperm*(txInd.f/g(phi)^2), c(2,3,1)),c(1,2),sum)
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
  omega <- drop(score1.vect(phiHat,lambdaHat) %*% (score3.vect(gammaHat) + p*l.vect(gammaHat))/n - 
                  sum(score3.vect(gammaHat) + p*l.vect(gammaHat))*apply(score1.vect(phiHat,lambdaHat),1,sum)/(n^2))   
  return(drop(solve(jack(phiHat,lambdaHat))[1:2,1:2] %*% omega)/(n*jack33))
}

# 'densRatio' computes maximum profile likelihood estimates of coefficients (and their variance estimates) in a mark density ratio model and returns a list containing:
#     'coef': estimates for alpha, beta, and lambda 
#     'var': the corresponding covariance matrix
#     'jack': the first two rows and columns of the limit estimating function in matrix form
#     'conv': a logical value indicating convergence of the estimating functions
# 'mark' is a numeric vector representing the mark variable, which is completely observed in all cases (i.e., failures)
# 'txInd' is the treatment group indicator (1 if treatment, 0 if control)
densRatio <- function(mark, txInd){
  V <- cbind(1,mark)
  z <- txInd
  nmark <- NCOL(V)
  ninf <- NROW(V)
  VV <- apply(V,1,tcrossprod)

  g <- function(theta){ exp(drop(V %*% theta)) }
  dG <- function(theta){ t(g(theta) * V) }
  d2G <- function(theta){ array(t(t(VV)*g(theta)),dim=c(nmark,nmark,ninf)) }
  dGdG <- function(theta){ array(apply(dG(theta),2,tcrossprod),dim=c(nmark,nmark,ninf)) }

  # profile score functions for the parameter of interest, theta, 
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

  return(list(coef=param.new, var=Var, jack=jack11(theta.new,lambda.new), conv=!(class(jackInv)=="try-error" | class(JackInv)=="try-error")))
}

#' Semiparametric Efficient Estimation and Testing for a Mark-Specific Proportional Hazards
#' Model with Multivariate Continuous Marks
#'
#' \code{sievePH} implements an efficient method of estimation for the multivariate mark-
#' specific hazard ratio in the competing risks failure time analysis framework in order
#' to assess mark-specific vaccine efficacy. The model is described in detail in Juraska
#' and Gilbert (2013) and improves efficiency of estimation by employing the semiparametric
#' method of maximum profile likelihood estimation in the vaccine-to-placebo mark density
#' ratio model. The model also enables the use of a more efficient estimation method,
#' proposed by Lu and Tsiatis (2008), for the overall log hazard ratio. The method employed
#' by \code{sievePH} is a complete-cases analysis of the mark in failures. 
#'
#' @param time a numeric vector specifying the observed time, defined as the minimum of the 
#' failure, censoring, and study time.
#' @param failInd a binary vector indicating the failure status (1 if failure, 0 if censored)
#' @param mark a numeric vector of the values of the mark variable, observed only in cases
#' @param txInd a binary vector indicating the treatment group (1 if treatment, 0 if control)
#'
#' @details
#' The conditional mark-specific hazard function can be factored into the product of the conditional 
#' mark density ratio and the ordinary marginal hazard function ignoring mark data. For the mark density 
#' ratio, following the assumption that the failure time and the mark variable are independent given the 
#' treatment group, a semiparametric density ratio model (Qin 1998) may be used. For the marginal hazard 
#' function, a Cox regression model is used.
#' 
#' Parameters in the mark density ratio are estimated with the maximum profile likelihood estimator, 
#' where the estimator is defined as the solution to the system of profile score functions for 
#' the parameter of interest and the Lagrange multiplier. The profile score functions are obtained by 
#' using the Lagrange multiplier method to maximize the semi-parametric log likelihood. This estimator 
#' is consistent and asymptotically normal. 
#' 
#' The parameter of interest in the marginal hazard function is estimated with the standard maximum 
#' partial likelihood estimator (MPLE) or the more efficient Lu and Tsiatis (2008) estimator, which 
#' leverages auxiliary data predictive of failure time (implemented in the R \code{speff2trial} package. 
#' Both estimators are consistent and asymptotically normal. 
#' 
#' The joint asymptotic distribution of the parameter estimators in the density ratio and Cox models 
#' is detailed in Juraska and Gilbert (2013) and is used to construct asymptotic pointwise Wald 
#' confidence intervals for the mark-specific vaccine efficacy.
#'
#' @return \code{sievePH} returns an object of class "sievePH" which can be processed by 
#' \code{\link{summary.sievePH}} to obtain or print a summary of the results. An object of class
#' "sievePH" is a list containing the following components:
#' \item{alphaHat}{the estimate for the \eqn{alpha} parameter in the density ratio model} 
#' \item{betaHat}{the estimate for the \eqn{beta} parameter in the density ratio model}
#' \item{gammaHat}{the estimate for the \eqn{gamma} parameter in the density ratio model}
#' \item{lambdaHat}{the estimate for \eqn{lambda}, the Lagrange multiplier utilized in 
#' the estimation procedure}
#' \item{cov}{the covariance matrix for \eqn{alpha}, \eqn{beta}, and \eqn{gamma}}
#' \item{ve}{a numeric vector of estimates of vaccine efficacy}
#' \item{mark}{a numeric vector of the values of the mark variable, observed only in cases}
#' \item{txInd}{a binary vector indicating the treatment group (1 if treatment, 0 if control)}
#' \item{nFail0}{the number of failures in the placebo group}
#' \item{nFail1}{the number of failures in the vaccine group}
#'
#' @examples
#'
#' @import survival
#'
#' @export
sievePH <- function(time, failInd, mark, txInd) {
  X <- time
  d <- failInd
  V <- mark
  Z <- txInd
  nFail0 <- sum(d*(1-Z))                             # number of failures in placebo group
  nFail1 <- sum(d*Z)                                 # number of failures in vaccine group

  dRatio <- densRatio(V[d==1],Z[d==1])

  if (dRatio$conv){

    # fit the Cox proportional hazards model to estimate the marginal hazard ratio
    phReg <- coxph(Surv(X,d)~Z)

    # parameter estimates
    thetaHat <- dRatio$coef
    gammaHat <- phReg$coef

    # variance and covariance estimates
    # order of columns in 'dRatio$var': alpha, beta1, beta2,...betak, lambda, where k is number of marks
    vthetaHat <- dRatio$var[1:2,1:2]
    vgammaHat <- drop(phReg$var)
    covThG <- covEst(X,d,V,Z,thetaHat[1:2],thetaHat[3],gammaHat)

    # vaccine efficacy estimate
    ve <- VE(V,thetaHat[1],thetaHat[2],gammaHat)

    # covariance matrix for alpha, beta1, gamma
    Sigma <- cbind(rbind(vthetaHat,covThG), c(covThG,vgammaHat))
    colnames(Sigma) <- rownames(Sigma) <- c("alpha", "beta1", "gamma")

    result <- list(mark = V, txInd = Z, nFail0 = nFail0, nFail1 = nFail1, alphaHat=thetaHat[1], betaHat=thetaHat[2], lambdaHat = thetaHat[3], gammaHat = gammaHat, 
                 ve = ve, cov = Sigma)
  } else {
    result <- list(mark = V, txInd = Z, nFail0 = nFail0, nFail1 = nFail1)
  }
  
  class(result) <- "sievePH"
  return(result)
}
