#' @import graphics
NULL
#' @import stats
NULL

# 'covEst' returns the estimated covariance matrix of 'phiHat' and 'lambdaHat' using Theorem 1 in Juraska and Gilbert (2013, Biometrics)
# 'eventTime' is the observed right-censored time on study
# 'find' is the failure indicator (0 if censored, 1 if failure)
# 'mark' is a data frame (with the same number of rows as the length of 'eventTime') specifying a multivariate mark (a numeric vector for a univariate mark is allowed), with NA for subjects with find=0.
# No missing mark values in subjects with find=1 are permitted.
# 'tx' is the treatment group indicator (1 if treatment, 0 if control)
# 'phiHat' is a vector of the alpha and beta estimates
# 'lambdaHat' is the estimate for lambda in the mark density ratio model
# 'gammaHat' is the estimate for gamma obtained in the marginal hazards model
covEst <- function(eventTime, find, mark, tx, phiHat, lambdaHat, gammaHat){
  # convert either a numeric vector or a data frame into a matrix
  mark <- as.matrix(mark)

  n <- length(eventTime)
  m <- sum(find)
  eventTime.f <- eventTime[find==1]
  V.f <- cbind(1,mark[find==1,])
  tx.f <- tx[find==1]
  eventTime.fM <- matrix(eventTime.f, nrow=n, ncol=m, byrow=TRUE)
  VV <- apply(V.f,1,tcrossprod)
  nmark <- NCOL(V.f)

  g <- function(phi){ exp(drop(V.f %*% phi)) }
  dG <- function(phi){ t(g(phi) * V.f) }
  d2G <- function(phi){ array(t(t(VV)*g(phi)),dim=c(nmark,nmark,m)) }
  dGdG <- function(phi){ array(apply(dG(phi),2,tcrossprod),dim=c(nmark,nmark,m)) }

  score1.vect <- function(phi, lambda){
    t((-lambda/(1+lambda*(g(phi)-1)) + tx.f/g(phi)) * t(dG(phi)))
  }
  xi <- function(gamma){ crossprod(eventTime>=eventTime.fM, tx*exp(gamma*tx)) }
  zeta <- function(gamma){ crossprod(eventTime>=eventTime.fM, exp(gamma*tx)) }
  eta <- drop(xi(gammaHat)/zeta(gammaHat))
  score3.vect <- function(gamma){ tx.f-eta }
  l.vect <- function(gamma){
    survprob.vect <- c(1, summary(survfit(Surv(eventTime,find)~1), times=sort(eventTime.f))$surv)
    surv.increm <- survprob.vect[-length(survprob.vect)] - survprob.vect[-1]
    eventTime.fMsq <- eventTime.fM[1:m,]
    crossprod(eventTime.f>=eventTime.fMsq, surv.increm*(tx.f*exp(gamma*tx.f) - eta*exp(gamma*tx.f))/zeta(gamma))
  }
  score1 <- function(phi, lambda){
    drop(-lambda * dG(phi) %*% (1/(1+lambda*(g(phi)-1))) + dG(phi) %*% (tx/g(phi)))
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
    term3 <- apply(aperm(d2Gperm*(tx.f/g(phi)), c(2,3,1)),c(1,2),sum)
    term4 <- apply(aperm(dGdGperm*(tx.f/g(phi)^2), c(2,3,1)),c(1,2),sum)
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

  omega <- drop(score1.vect(phiHat,lambdaHat) %*% (score3.vect(gammaHat) + p*l.vect(gammaHat))/n -
                  sum(score3.vect(gammaHat) + p*l.vect(gammaHat))*apply(score1.vect(phiHat,lambdaHat),1,sum)/(n^2))
  return(drop(solve(jack(phiHat,lambdaHat))[1:nmark,1:nmark] %*% omega)/(n*jack33))
}

# 'densRatio' computes maximum profile likelihood estimates of coefficients (and their variance estimates) in a mark density ratio model and returns a list containing:
#     'coef': estimates for alpha, beta, and lambda
#     'var': the corresponding covariance matrix
#     'jack': the first two rows and columns of the limit estimating function in matrix form
#     'conv': a logical value indicating convergence of the estimating functions
# 'mark' is a data frame representing a multivariate mark variable (a numeric vector for a univariate mark is allowed), which is completely observed in all cases (i.e., failures). No missing mark values are permitted.
# 'tx' is the treatment group indicator (1 if treatment, 0 if control)
densRatio <- function(mark, tx){
  # convert either a numeric vector or a data frame into a matrix
  mark <- as.matrix(mark)

  V <- cbind(1,mark)
  z <- tx
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
    if (!inherits(jackInv, "try-error")){
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
  if (!inherits(JackInv, "try-error")){
    Var <- ninf * JackInv %*% SigmaHat(theta.new,lambda.new) %*% JackInv
    names(param.new) <- rownames(Var) <- colnames(Var) <- c("alpha",
                                                            paste("beta",1:(nmark-1),sep=""),"lambda")
  } else {
    Var <- NULL
  }

  return(list(coef=param.new, var=Var, jack=jack11(theta.new,lambda.new), conv=!(inherits(jackInv, "try-error") | inherits(JackInv, "try-error"))))
}

#' Semiparametric Estimation of Coefficients in a Mark-Specific Proportional Hazards
#' Model with a Multivariate Continuous Mark, Fully Observed in All Failures
#'
#' \code{sievePH} implements the semiparametric estimation method of Juraska and Gilbert (2013) for the multivariate mark-
#' specific hazard ratio in the competing risks failure time analysis framework. It employs (i) the semiparametric
#' method of maximum profile likelihood estimation in the treatment-to-placebo mark density
#' ratio model (Qin, 1998) and (ii) the ordinary method of maximum partial likelihood estimation of the overall log hazard ratio in the Cox model.
#' \code{sievePH} requires that the multivariate mark data are fully observed in all failures.
#'
#' @param eventTime a numeric vector specifying the observed right-censored time to the event of interest
#' @param eventInd a numeric vector indicating the event of interest (1 if event, 0 if right-censored)
#' @param mark either a numeric vector specifying a univariate continuous mark or a data frame specifying a multivariate continuous mark.
#' No missing values are permitted for subjects with \code{eventInd = 1}. For subjects with \code{eventInd = 0}, the value(s) in \code{mark} should be set to \code{NA}.
#' @param tx a numeric vector indicating the treatment group (1 if treatment, 0 if placebo)
#'
#' @details
#' \code{sievePH} considers data from a randomized placebo-controlled treatment efficacy trial with a time-to-event endpoint.
#' The parameter of interest, the mark-specific hazard ratio, is the ratio (treatment/placebo) of the conditional mark-specific hazard functions.
#' It factors as the product of the mark density ratio (treatment/placebo) and the ordinary marginal hazard function ignoring mark data.
#' The mark density ratio is estimated using the method of Qin (1998), while the marginal hazard ratio is estimated using \code{coxph()} in the \code{survival} package.
#' Both estimators are consistent and asymptotically normal. The joint asymptotic distribution of the estimators is detailed in Juraska and Gilbert (2013).
#'
#' @return An object of class \code{sievePH} which can be processed by
#' \code{\link{summary.sievePH}} to obtain or print a summary of the results. An object of class
#' \code{sievePH} is a list containing the following components:
#' \itemize{
#' \item \code{DRcoef}: a numeric vector of estimates of coefficients \eqn{\phi} in the weight function \eqn{g(v, \phi)} in the density ratio model
#' \item \code{DRlambda}: an estimate of the Lagrange multiplier in the profile score functions for \eqn{\phi} (that arises by profiling out the nuisance parameter)
#' \item \code{DRconverged}: a logical value indicating whether the estimation procedure in the density ratio model converged
#' \item \code{logHR}: an estimate of the marginal log hazard ratio from \code{coxph()} in the \code{survival} package
#' \item \code{cov}: the estimated joint covariance matrix of \code{DRcoef} and \code{logHR}
#' \item \code{coxphFit}: an object returned by the call of \code{coxph()}
#' \item \code{nPlaEvents}: the number of events observed in the placebo group
#' \item \code{nTxEvents}: the number of events observed in the treatment group
#' \item \code{mark}: the input object
#' \item \code{tx}: the input object
#' }
#'
#' @references Juraska, M. and Gilbert, P. B. (2013), Mark-specific hazard ratio model with multivariate continuous marks: an application to vaccine efficacy. \emph{Biometrics} 69(2):328–337.
#'
#' Qin, J. (1998), Inferences for case-control and semiparametric two-sample density ratio models. \emph{Biometrika} 85, 619–630.
#'
#' @examples
#' n <- 500
#' tx <- rep(0:1, each=n/2)
#' tm <- c(rexp(n/2, 0.2), rexp(n/2, 0.2 * exp(-0.4)))
#' cens <- runif(n, 0, 15)
#' eventTime <- pmin(tm, cens, 3)
#' eventInd <- as.numeric(tm <= pmin(cens, 3))
#' mark1 <- ifelse(eventInd==1, c(rbeta(n/2, 2, 5), rbeta(n/2, 2, 2)), NA)
#' mark2 <- ifelse(eventInd==1, c(rbeta(n/2, 1, 3), rbeta(n/2, 5, 1)), NA)
#'
#' # fit a model with a univariate mark
#' fit <- sievePH(eventTime, eventInd, mark1, tx)
#'
#' # fit a model with a bivariate mark
#' fit <- sievePH(eventTime, eventInd, data.frame(mark1, mark2), tx)
#'
#' @seealso \code{\link{summary.sievePH}}, \code{\link{plot.summary.sievePH}}, \code{\link{testIndepTimeMark}} and \code{\link{testDensRatioGOF}}
#'
#' @import survival
#'
#' @export
sievePH <- function(eventTime, eventInd, mark, tx){
  if (is.numeric(mark)){ mark <- data.frame(mark) }

  nPlaEvents <- sum(eventInd * (1-tx))
  nTxEvents <- sum(eventInd * tx)

  dRatio <- densRatio(mark[eventInd==1, ], tx[eventInd==1])

  # fit the Cox proportional hazards model to estimate the marginal hazard ratio
  phReg <- survival::coxph(Surv(eventTime, eventInd) ~ tx)

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
    covThG <- covEst(eventTime, eventInd, mark, tx, thetaHat[-lastComp], thetaHat[lastComp], gammaHat)

    # covariance matrix for alpha, beta1, beta2,..., betak, gamma
    Sigma <- cbind(rbind(vthetaHat, covThG), c(covThG, vgammaHat))
    colnames(Sigma) <- rownames(Sigma) <- c("alpha", paste0("beta", 1:NCOL(mark)), "gamma")

    out$DRcoef <- thetaHat[-lastComp]
    out$DRlambda <- thetaHat[lastComp]
    out$cov <- Sigma
  } else {
    warning("The estimation method in the density ratio model did not converge.")
  }

  class(out) <- "sievePH"
  return(out)
}
