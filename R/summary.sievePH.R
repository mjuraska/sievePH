# 'profileLRtest' performs the profile likelihood ratio test and returns a list consisting of the test statistics and p-values
# 'mark' is a data frame specifying a multivariate mark (a numeric vector for a univariate mark is allowed), which is completely observed in all cases (i.e., failures). No missing mark values are permitted.
# 'tx' is a numeric vector indicating the treatment group (1 if treatment, 0 if control)
profileLRtest <- function(mark, tx, thetaHat, lambdaHat){
  mark <- as.matrix(mark)
  nMark <- NCOL(mark)

  g <- function(theta, mark){ exp(drop(cbind(1, mark) %*% theta)) }
  loglik <- function(theta, lambda, mark){ -sum(log(1 + lambda * (g(theta, mark) - 1))) + sum(tx * log(g(theta, mark))) }

  teststat <- 2 * (loglik(thetaHat, lambdaHat, mark) - loglik(rep(0, nMark + 1), 0, mark))
  pval <- 1 - pchisq(teststat, nMark)
  return(list(teststat=teststat, pval=pval))
}

### univariate two-sided Wald test of H0: 'est'=0, where 'est' is alpha, gamma, or a component of beta
waldH0.2sided.pval <- function(est, vEst){
  testStat <- est / sqrt(vEst)
  return(2*pnorm(-abs(testStat)))
}

#' Summarizing Mark-Specific Proportional Hazards Model Fits
#'
#' \code{summary} method for an object of class \code{sievePH}.
#'
#' @aliases print.summary.sievePH
#' @param object an object of class \code{sievePH}, usually a result of a call to \code{\link{sievePH}}
#' @param markGrid a matrix specifying a grid of multivariate mark values, where rows correspond to different values on the (multivariate) grid and columns correspond to components of the mark. A numeric vector is allowed
#' for univariate marks. The point and interval estimates of the \code{contrast} are calculated on this grid.
#' @param contrast a character string specifying the treatment effect parameter of interest. The default value is \code{"te"} (treatment efficacy); other options are \code{"hr"} (hazard ratio) and \code{"loghr"} (log hazard ratio).
#' @param sieveAlternative a character string specifying the alternative hypothesis for the sieve tests, which can be either \code{"twoSided"} (default) or, in case of a univariate mark, \code{"HR.decrease"} or \code{"HR.increase"}.
#' The one-sided option is unavailable for a multivariate mark.
#' @param confLevel the confidence level (0.95 by default) of reported confidence intervals
#'
#' @details
#' \code{print.summary.sievePH} prints a formatted summary of results. Inference about coefficients in the mark-specific proportional hazards model is tabulated. Additionally, a summary is generated
#' from the likelihood-ratio and Wald tests of two relevant null hypotheses: (1) \{\eqn{H_0: HR(v)=1} for all \eqn{v}\}, and (2) \{\eqn{H_0: HR(v)=HR} for all \eqn{v}\}. For the tests of (2) and a univariate
#' mark, \code{sieveAlternative} controls the choice of the alternative hypothesis.
#'
#' @return
#' An object of class \code{summary.sievePH}, which is a list with the following components:
#' \itemize{
#' \item \code{coef}: a data frame summarizing point and interval estimates of the density ratio model coefficients and the marginal log hazard ratio (the confidence level is specified by \code{confLevel}), and p-values from the
#' two-sided Wald test of the null hypothesis that the parameter equals zero
#' \item \code{pLR.HRunity.2sided}: a numeric vector with two named components: \code{pLR.dRatio.2sided} is a p-value from the two-sided profile likelihood-ratio test of the null hypothesis \eqn{H_0: \beta=0}, where \eqn{\beta} is the
#' vector of mark coefficients in the mark density ratio model, and \code{pLR.cox.2sided} is a p-value from the two-sided partial likelihood-ratio test of the null hypothesis \eqn{H_0: \gamma=0}, where \eqn{\gamma} is the
#' marginal log hazard ratio in the Cox model. The two p-values are intended for the use of the Simes (1986) procedure as described on page 4 in Juraska and Gilbert (2013).
#' \item \code{pWald.HRunity.2sided}: a p-value from the two-sided Wald test of the null hypothesis \{\eqn{H_0: HR(v)=1} for all \eqn{v}\}
#' \item \code{pWtWald.HRunity.1sided}: a p-value from the one-sided weighted Wald test of the null hypothesis \{\eqn{H_0: HR(v)=1} for all \eqn{v}\} against the alternative hypothesis \{\eqn{H_1: HR < 1} and \eqn{HR(v)} is
#' increasing in each component of \eqn{v}\}
#' \item \code{pLR.HRconstant.2sided}: a p-value from the two-sided profile likelihood-ratio test of the null hypothesis \{\eqn{H_0: HR(v)=HR} for all \eqn{v}\}. This component is available if \code{sieveAlternative="twoSided"}.
#' \item \code{pLR.HRconstant.1sided}: a numeric vector with two named components: \code{pLR.dRatio.2sided} is a p-value from the two-sided profile likelihood-ratio test of the null hypothesis \{\eqn{H_0: HR(v)=HR} for all \eqn{v}\},
#' and \code{estBeta} is the point estimate of the univariate mark coefficient in the density ratio model. This component is available if the mark is univariate and \code{sieveAlternative="oneSided"}.
#' \item \code{pWald.HRconstant.2sided}: a p-value from the two-sided Wald test of the null hypothesis \{\eqn{H_0: HR(v)=HR} for all \eqn{v}\}. This component is available if \code{sieveAlternative="twoSided"}.
#' \item \code{pWald.HRconstant.1sided}: a p-value from the one-sided Wald test of the null hypothesis \{\eqn{H_0: HR(v)=HR} for all \eqn{v}\} against the alternative hypothesis \{\eqn{H_1: HR(v)} is increasing in \eqn{v}\} or \{\eqn{H_1: HR(v)} is decreasing in \eqn{v}\}.
#' This component is available if the mark is univariate and \code{sieveAlternative="HRdecrease"} or \code{sieveAlternative="HRincrease"}.
#' \item \code{te}: a data frame summarizing point and interval estimates of the mark-specific treatment efficacy on the grid of mark values in \code{markGrid} (available if \code{contrast="te"}). The confidence level is specified
#' by \code{confLevel}.
#' \item \code{hr}: a data frame summarizing point and interval estimates of the mark-specific hazard ratio on the grid of mark values in \code{markGrid} (available if \code{contrast="hr"}). The confidence level is specified by
#' \code{confLevel}.
#' \item \code{loghr}: a data frame summarizing point and interval estimates of the mark-specific log hazard ratio on the grid of mark values in \code{markGrid} (available if \code{contrast="loghr"}). The confidence level is specified by
#' \code{confLevel}.
#' }
#'
#' @references Juraska, M. and Gilbert, P. B. (2013), Mark-specific hazard ratio model with multivariate continuous marks: an application to vaccine efficacy. \emph{Biometrics} 69(2):328–337.
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
#' # fit a model with a bivariate mark
#' fit <- sievePH(eventTime, eventInd, data.frame(mark1, mark2), tx)
#' sfit <- summary(fit, markGrid=matrix(c(0.3, 0.3, 0.6, 0.3, 0.3, 0.6, 0.6, 0.6),
#'                                      ncol=2, byrow=TRUE))
#' # print the formatted summary
#' sfit
#' # treatment efficacy estimates on the grid
#' sfit$te
#'
#' @seealso \code{\link{sievePH}}
#'
#' @export
summary.sievePH <- function(object, markGrid,
                            contrast = c("te", "hr", "loghr"),
                            sieveAlternative = c("twoSided","HRincrease", "HRdecrease"), confLevel = 0.95, ...){

  if (missing(markGrid)){ stop("The grid of mark values in 'markGrid' is missing.") }

  contrast <- match.arg(contrast, choices = c("te", "hr", "loghr"))
  sieveAlternative <- match.arg(sieveAlternative, choices = c("twoSided","HRincrease", "HRdecrease"))

  nMark <- NCOL(object$mark)
  alphaHat <- object$DRcoef[1]
  betaHat <- object$DRcoef[-1]
  lambdaHat <- object$DRlambda
  gammaHat <- object$logHR
  variances <- diag(object$cov)
  vAlphaHat <- variances[1]
  vBetaHat <- variances[2:(nMark + 1)]
  vGammaHat <- variances[nMark + 2]

  phRegSummary <- summary(object$coxphFit)

  # probability level to be used for quantiles in the construction of confidence bounds
  probLevel <- 1 - (1-confLevel)/2

  # 'sievePH' requires that the mark data be fully observed in all failures
  isNA <- attr(na.omit(object$mark), "na.action")

  ### two-sided profile likelihood ratio test of constancy of HR (H0: beta=0)
  pLR.dRatio.2sided <- profileLRtest(object$mark[-isNA, ], object$tx[-isNA], c(alphaHat, betaHat), lambdaHat)$pval

  ### two-sided likelihood ratio test of H0: HR(v)=1 (i.e., HR unity) against H1: HR(v)!=1 intended for the use of the Simes (1986) procedure as described on page 4 in Juraska and Gilbert (2013, Biometrics)
  ### H0 is equivalent to H0*: beta=0 and gamma=0, and H1 is equivalent to H1*: beta!=0 or gamma!=0
  ### a named vector with two two-sided p-values, one from the profile LR test for beta and one from the partial LR test for gamma
  ### the components of the vector are named 'pLR.dRatio.2sided' and 'pLR.cox.2sided'
  pLR.HRunity.2sided <- c(pLR.dRatio.2sided, pLR.cox.2sided = phRegSummary$logtest[3])
  names(pLR.HRunity.2sided) <- c("pLR.dRatio.2sided", "pLR.cox.2sided")

  ### two-sided Wald test of H0: HR(v)=1 against H1: HR(v)!=1
  ### H0 is equivalent to H0*: beta=0 and gamma=0, and H1 is equivalent to H1*: beta!=0 or gamma!=0
  pWald.HRunity.2sided <- 1 - pchisq(drop(t(c(betaHat, gammaHat)) %*% solve(object$cov[-1, -1]) %*% c(betaHat, gammaHat)), nMark + 1)

  ### one-sided weighted Wald-type test of H0: HR(v)=1 against H1: HR(v)<1 and HR(v) increasing in each component of v
  ### the test statistic is in formula (11) on page 4 in Juraska and Gilbert (2013, Biometrics)
  if (nMark==1){
    weighted.waldH00 <- (betaHat / vBetaHat - gammaHat / vGammaHat) / sqrt(1 / vBetaHat + 1 / vGammaHat - 2 * object$cov[2, 3] / (vBetaHat * vGammaHat))
  } else {
    # 'covBeta' is always at leat a 2x2 matrix
    covBeta <- object$cov[2:(nMark + 1), 2:(nMark + 1)]
    covBeta[lower.tri(covBeta, diag=TRUE)] <- 0
    weighted.waldH00 <- (sum(betaHat / vBetaHat) - gammaHat / vGammaHat) / sqrt(sum(1 / vBetaHat) + 1 / vGammaHat + 2 * drop(t(1 / vBetaHat) %*% covBeta %*% (1 / vBetaHat)) -
                                                                                  2 * (1 / vGammaHat) * drop(t(1 / vBetaHat) %*% object$cov[2:(nMark + 1), nMark + 2]))
  }
  pWtWald.HRunity.1sided <- 1 - pnorm(weighted.waldH00)

  if (sieveAlternative %in% c("HRdecrease", "HRincrease") & nMark > 1){ warning("One-sided sieve tests are available for univariate marks only.") }

  if (sieveAlternative=="HRdecrease" & nMark==1){

    ### 1-sided Wald test of H0: HR(v)=HR (i.e., beta=0) vs H1: HR(v) decreasing in v (i.e., beta<0)
    waldH0 <- betaHat / sqrt(vBetaHat)
    pWald.HRconstant <- pnorm(waldH0)
    
    # ### The labels are 'pWald.beta.2sided' and 'estBeta1', 'estBeta2', etc. (if the dimension of beta is 1, then only 'estBeta')
    # pWald.HRconstant <- 1 - pchisq(drop(t(betaHat) %*% solve(object$cov[-c(1, length(variances)), -c(1, length(variances))]) %*% betaHat), nMark)
    # names(pWald.HRconstant) <- c("pWald.beta.2sided", ifelse(nMark==1, "estBeta", sapply(1:nMark, function(x) paste0("estBeta", x))))

    ### 1-sided likelihood ratio test of H0: HR(v)=HR (i.e., beta=0) vs H1: HR(v) increasing in v (i.e., beta>0)
    ### A named vector with the following components: the two-sided profile LR test p-value, and the point estimate of beta
    ### the labels are 'pLR.dRatio.2sided' and 'estBeta'
    pLR.HRconstant <- c(pLR.dRatio.2sided, betaHat)
    names(pLR.HRconstant) <- c("pLR.dRatio.2sided", "estBeta")

  }else if(sieveAlternative=="HRincrease" & nMark==1){
    ### 1-sided Wald test of H0: HR(v)=HR (i.e., beta=0) vs H1: HR(v) increasing in v (i.e., beta>0)
    waldH0 <- betaHat / sqrt(vBetaHat)
    pWald.HRconstant <- 1 - pnorm(waldH0)
    
    pLR.HRconstant <- c(pLR.dRatio.2sided, betaHat)
    names(pLR.HRconstant) <- c("pLR.dRatio.2sided", "estBeta")
  }else {

    ### 2-sided Wald test of H0: HR(v)=HR (i.e., beta=0) vs H1: HR(v)!=HR (i.e., beta!=0)
    pWald.HRconstant <- 1 - pchisq(drop(t(betaHat) %*% solve(object$cov[2:(nMark + 1), 2:(nMark + 1)]) %*% betaHat), nMark)

    ### 2-sided likelihood ratio test of H0: HR(v)=HR (i.e., beta=0) vs H1: HR(v)!=HR (i.e., beta!=0)
    pLR.HRconstant <- pLR.dRatio.2sided
  }

  # initialize component 'coef' of the output list
  coef <- as.data.frame(matrix(nrow=nMark + 2, ncol=4))
  colnames(coef) <- c("Estimate", "LB", "UB", "pWald")
  rownames(coef) <- c("DR Intercept", colnames(object$mark), "Marginal Log HR")

  est <- c(alphaHat, betaHat, gammaHat)
  vEst <- c(vAlphaHat, vBetaHat, vGammaHat)

  coef[, "Estimate"] <- est
  coef[, c("LB", "UB")] <- est + qnorm(probLevel) * sqrt(vEst) %o% c(-1, 1)
  coef[, "pWald"] <- c(waldH0.2sided.pval(est[-length(est)], vEst[-length(vEst)]), phRegSummary$waldtest[3])

  # linear score estimates and variance estimates
  lScore <- drop(cbind(1, markGrid, 1) %*% c(alphaHat, betaHat, gammaHat))
  varLscore <- diag(cbind(1, markGrid, 1) %*% object$cov %*% t(cbind(1, markGrid, 1)))
  # the below 'varLscore' for univariate marks is identical to the above 'varLscore'
  # varLscore <- vAlphaHat + (markGrid^2) * vBetaHat + vGammaHat + 2 * markGrid * object$cov[1, 2] + 2 * object$cov[1, 3] + 2 * markGrid * object$cov[2, 3]

  # linear score lower and upper confidence limits with confidence level specified by 'confLevel'
  lLB <- lScore - sqrt(varLscore) * qnorm(probLevel)
  lUB <- lScore + sqrt(varLscore) * qnorm(probLevel)

  # contrast estimates and confidence bounds
  if (contrast=="te"){
    est <- 1 - exp(lScore)
    lb <- 1 - exp(lUB)
    ub <- 1 - exp(lLB)
  } else if (contrast=="hr"){
    est <- exp(lScore)
    lb <- exp(lLB)
    ub <- exp(lUB)
  } else if (contrast=="loghr"){
    est <- lScore
    lb <- lLB
    ub <- lUB
  }

  # assemble the contrast matrix of the output list
  contrastDF <- data.frame(markGrid, est, lb, ub)
  colnames(contrastDF) <- c(colnames(object$mark), switch(contrast, te="TE", hr="HR", loghr="LogHR"), "LB", "UB")
  out <- list(coef, pLR.HRunity.2sided, pWald.HRunity.2sided, pWtWald.HRunity.1sided, pLR.HRconstant, 
              pWald.HRconstant, contrastDF)
  names(out) <- c("coef", "pLR.HRunity.2sided", "pWald.HRunity.2sided", "pWtWald.HRunity.1sided",
                  paste0("pLR.HRconstant.", ifelse(sieveAlternative %in% c("HRdecrease", "HRincrease") & nMark==1, "1", "2"), "sided"),
                  paste0("pWald.HRconstant.", ifelse(sieveAlternative %in% c("HRdecrease", "HRincrease") & nMark==1, "1", "2"), "sided"),
                  contrast)
  

  class(out) <- "summary.sievePH"
  return(out)
}

#' @rdname summary.sievePH
#' @param x an object of class \code{summary.sievePH}, usually a result of a call to \code{summary.sievePH}
#' @param digits the number of significant digits to use when printing (4 by default)
#' @param ... further arguments passed to or from other methods
#' @export
print.summary.sievePH <- function(x, digits=4, ...){
  cat("\nCoefficients:\n")
  print(x$coef, digits=digits, print.gap=2)
  cat("\n")
  cat("Tests of H0: HR(v) = 1 for all v:\n")
  cat("Two-sided likelihood-ratio test:\n")
  cat("  Density ratio model profile likelihood-ratio test p-value: ", format(x$pLR.HRunity.2sided["pLR.dRatio.2sided"], digits=digits, nsmall=digits), "\n", sep="")
  cat("  Cox model partial likelihood-ratio test p-value: ", format(x$pLR.HRunity.2sided["pLR.cox.2sided"], digits=digits, nsmall=digits), "\n", sep="")
  cat("Two-sided Wald test p-value: ", format(x$pWald.HRunity.2sided, digits=digits, nsmall=digits), "\n", sep="")
  cat("One-sided weighted Wald test p-value: ", format(x$pWtWald.HRunity.1sided, digits=digits, nsmall=digits), "\n\n", sep="")
  cat("Tests of H0: HR(v) = HR for all v:\n")
  if (is.null(x$pLR.HRconstant.2sided)){
    cat("Two-sided likelihood-ratio test p-value: ", format(x$pLR.HRconstant.1sided["pLR.dRatio.2sided"], digits=digits, nsmall=digits), "\n", sep="")
    cat("  Point estimate of the mark coefficient: ", format(x$pLR.HRconstant.1sided["estBeta"], digits=digits, nsmall=digits), "\n", sep="")
    cat("One-sided Wald test p-value: ", format(x$pWald.HRconstant.1sided, digits=digits, nsmall=digits), "\n", sep="")
  } else {
    cat("Two-sided likelihood-ratio test p-value: ", format(x$pLR.HRconstant.2sided, digits=digits, nsmall=digits), "\n", sep="")
    cat("Two-sided Wald test p-value: ", format(x$pWald.HRconstant.2sided, digits=digits, nsmall=digits), "\n", sep="")
  }
}


