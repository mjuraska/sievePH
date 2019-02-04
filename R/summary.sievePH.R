# 'VE' returns vaccine efficacy values given parameters for the grid of mark values, V,
# and values for alpha, beta, and gamma
VE <- function(v, alpha, beta, gamma){ 1 - exp(alpha + beta*v + gamma) }

dVE <- function(v, alpha, beta, gamma){ -exp(alpha+beta*v+gamma)*c(1,v,1) }

# 'seVE' calculates the standard error given parameters for the grid of mark values, v, values for alpha, beta, and gamma, and the corresponding covariance matrix
seVE <- function(v, cov, alpha, beta, gamma){
  sapply(v, function(mark){ drop(sqrt(t(dVE(mark,alpha,beta,gamma)) %*% cov %*% dVE(mark,alpha,beta,gamma))) })
}

# 'profileLRtest' performs the profile likelihood ratio test and returns a list consisting of the test statistics and p-values
# 'mark' is a matrix specifying a multivariate mark (a numeric vector for a univariate mark is allowed), which is completely observed in all cases (i.e., failures). No missing mark values are permitted.
# 'tx' is a binary vector indicating the treatment group (1 if treatment, 0 if control)
profileLRtest <- function(mark, tx, thetaHat, lambdaHat){
  mark <- as.matrix(mark)
  nMark <- NCOL(mark)

  g <- function(theta, mark){ exp(drop(cbind(1, mark) %*% theta)) }
  loglik <- function(theta, lambda, mark){ -sum(log(1 + lambda * (g(theta, mark) - 1))) + sum(tx * log(g(theta, mark))) }
  
  teststat <- 2 * (loglik(thetaHat, lambdaHat, mark) - loglik(rep(0, nMark + 1), 0, mark))
  pval <- 1 - pchisq(teststat, nMark)
  return(list(teststat=teststat, pval=pval))
}

#' Summarizing Results for Semiparametric Efficient Estimation and Testing for a Multivariate Mark-Specific Hazard Ratio Model
#'
#' summary method for an object of class "sievePH"
#'
#' @param object an object of class "sievePH", usually a result of a call to \code{\link{sievePH}}.
#' @param markGrid a matrix specifying a grid of multivariate mark values, where rows correspond to different values on the (multivariate) grid and columns correspond to components of the mark. The point and interval
#' estimates of the \code{contrast} are calculated on the grid.
#' @param contrast a string specifying the estimate of interest. Default is \code{"ve"}; other options are \code{"hr"} and \code{"loghr"}.
#' @param sieveAlternative a string specifying the type of alternative for the sieve test, with possible values being \code{"twoSided"} (default) or \code{"oneSided").
#' @param confLevel the confidence level to be used for the reported confidence intervals. Default is \code{0.95}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' \code{print.summary.speff} prints a formatted summary of results. An inferential table is
#' produced with point and interval estimates for the mark-specific hazard ratio, standard
#' error estimates, and Wald test p-values for the two relevant hypothesis tests of no
#' vaccine efficacy and no sieve effect.
#'
#' The user can specify whether a one-sided or two-sided hypothesis test is to be performed
#' when determining if vaccine efficacy varies by viral divergence of the mark
#' (i.e., null hypothesis H0:HR(v) = HR).
#'
#' @return
#' A list with the following components:
#' \itemize{
#' \item \code{coef}: a numeric matrix
#' \item \code{pLR.HRunity.2sided}: the p-value for the two-sided likelihood ratio test of unity
#' \item \code{pWald.HRunity.2sided}: the p-value for the two-sided WAld test of unity
#' \item \code{pWtWald.HRunity.1sided}: the p-value for the one-sided weighted Wald test of unity
#' \item \code{pLR.HRconstant.2sided}: the p-value for the two-sided likelihood ratio test of constancy
#' \item \code{pWald.HRconstant.2sided}: the p-value for the two-sided Wald test of constancy
#' \item \code{ve}: a numeric matrix containing the values of the multivariate mark variable and point and interval estimates for the component specified by the \code{contrast} input parameter, with the confidence level specified by the \code{confLevel} input parameter
#' }
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
#' # fit a model with a bivariate mark
#' fit <- sievePH(eventTime, eventInd, data.frame(mark1, mark2), tx)
#' summary(fit)
#'
#' @seealso \code{\link{sievePH}}
#'
#' @export
summary.sievePH <- function(object, markGrid,
                            contrast = c("ve", "hr", "loghr"),
                            sieveAlternative = c("twoSided","oneSided"), confLevel = 0.95, ...){

  if (missing(markGrid)){ stop("The grid of mark values in 'markGrid' is missing.") }

  sieveAlternative <- match.arg(sieveAlternative, choices = c("twoSided","oneSided"))
  contrast <- match.arg(contrast, choices = c("ve", "hr", "loghr"))

  nMark <- NCOL(object$mark)
  betaHat <- object$DRcoef[-1]
  alphaHat <- object$DRcoef[1]
  lambdaHat <- object$DRlambda
  gammaHat <- object$logHR
  variances <- diag(object$cov)
  vAlphaHat <- variances[1]
  vBetaHat <- variances[2:(nMark+1)]
  vGammaHat <- variances[length(variances)]

  phRegSummary <- summary(object$coxphFit)

  # quantile to be used in confidence bounds
  quantile <- 1 - (1-confLevel)/2

  # 'sievePH' requires that the mark data be fully observed in all failures
  isNA <- attr(na.omit(object$mark), "na.action")

  ### two-sided profile likelihood ratio test of constancy of HR (H0: beta=0)
  pLR.dRatio.2sided <- profileLRtest(object$mark[-isNA, ], object$tx[-isNA], c(alphaHat, betaHat), lambdaHat)$pval

  ### two-sided likelihood ratio test of H0: HR(v)=1 (i.e., HR unity) against H1: HR(v)!=1 intended for the use of the Simes (1986) procedure as described on page 4 in Juraska and Gilbert (2013, Biometrics)
  ### H0 is equivalent to H0*: beta=0 and gamma=0, and H1 is equivalent to H1*: beta!=0 or gamma!=0
  ### a named vector with two two-sided p-values, one from the profile LR test for beta and one from the partial LR test for gamma
  ### the components of the vector are named 'pLR.dRatio.2sided' and 'pLR.cox.2sided'
  pLR.HRunity.2sided <- c(pLR.dRatio.2sided, pLR.cox.2sided = phRegSummary$logtest[3])

  ### two-sided Wald test of H0: HR(v)=1 against H1: HR(v)!=1
  ### H0 is equivalent to H0*: beta=0 and gamma=0, and H1 is equivalent to H1*: beta!=0 or gamma!=0
  pWald.HRunity.2sided <- 1 - pchisq(drop(t(c(betaHat, gammaHat)) %*% solve(object$cov[-1, -1]) %*% c(betaHat, gammaHat)), nMark + 1)

  ### one-sided weighted Wald-type test of H0: HR(v)=1 against H1: HR(v)<1 and HR(v) increasing in each component of v
  weighted.waldH00 <- (betaHat/vBetaHat - gammaHat/vGammaHat) / sqrt(1/vBetaHat + 1/vGammaHat - 2*object$cov[3,2]/(vBetaHat*vGammaHat))
  pWtWald.HRunity.1sided <- 1 - pnorm(weighted.waldH00)

  ### univariate two-sided Wald test of H0: 'est'=0, where 'est' is alpha, gamma, or a component of beta
  waldH0.2sided.pval <- function(est, vEst) {
    testStat <- est / sqrt(vEst)
    return(2*pnorm(-abs(testStat)))
  }

  if (sieveAlternative == "oneSided") {

    ### 1-sided Wald test of H0: HR(v)=HR (beta=0) vs alternative that beta > 0
    waldH0 <- betaHat/sqrt(vBetaHat)  ### still univariate
    pWald.HRconstant <- 1 - pnorm(waldH0)

    ### 1-sided likelihood ratio test of H0: HR(v)=HR (beta=0) vs alternative that beta > 0
    ### A named vector with the following components: the two-sided profile LR test p-value, and point estimates of the components of the vector beta
    ### the labels are 'pLR.beta.2sided' and 'estBeta1', 'estBeta2', etc. (if the dimension of beta is 1, then only 'estBeta')
    pLR.HRconstant <- c(pLR.beta.2sided, betaHat)
    names(pLR.HRconstant) <- c("pLR.beta.2sided", ifelse(nMark==1, "estBeta", sapply(1:nMark, function(x) paste0("estBeta", x))))

  } else {

    ### 2-sided Wald test of H0: HR(v)=HR (beta=0)
    pWald.HRconstant <- 1 - pchisq(drop(t(betaHat) %*% solve(object$cov[-c(1, length(variances)), -c(1, length(variances))]) %*% betaHat), nMark)
    
    ### 2-sided likelihood ratio test of H0: HR(v)=HR (beta=0)
    pLR.HRconstant <- pLR.beta.2sided
  }


  # create output matrix of coefficients
  coef <- matrix(nrow = nMark + 2, ncol = 5)
  colnames(coef) <- c("Estimate", "LB", "UB", "pLR", "pWald")
  rownames(coef) <- c("DR Intercept", colnames(object$mark), "Marginal Log HR")

  est <- c(alphaHat, betaHat, gammaHat)
  vEst <- c(vAlphaHat, vBetaHat, vGammaHat)

  coef[, "Estimate"] <- est
  coef[, c("LB", "UB")] <- t(mapply(function(x, y) {x + qnorm(quantile)*sqrt(y) %o% c(-1,1)}, est, vEst))
  coef[, "pLR"] <- c(rep(pLR.beta.2sided, nMark + 1), phRegSummary$logtest[3])  #### change?
  coef[, "pWald"] <- c(waldH0.2sided.pval(est[-length(est)], vEst[-length(vEst)]), phRegSummary$waldtest[3])

  # linear score estimates and variance estimates
  lScore <- log(1-ve)
  if (nMark > 1) {
    combs <- combn(nMark, 2)  # distinct combinations of indices for betaHat
    covAlphaBeta <- sapply(1:nMark, function(i) {2*markGrid[,i]*object$cov[1, i+1]})  # covariance of alpha and elements of betaHat
    covBeta <- sapply(1:ncol(combs), function(i) {2*markGrid[,combs[1,i]]*markGrid[,combs[2,i]]*object$cov[combs[1,i]+1, combs[2,i]+1]})  # covariance of elements of betaHat
    varLscore <- vAlphaHat + t(t(markGrid^2)*vBetaHat) + vGammaHat + covAlphaBeta + covBeta
  } else {  # univariate mark
    varLscore <- vAlphaHat + (markGrid^2)*vBetaHat + 2*markGrid*object$cov[1,2] + vGammaHat
  }
  # linear score lower and upper confidence limits with confidence level specified by 'confLevel'
  lLB <- lScore - (sqrt(varLscore)*qnorm(quantile))
  lUB <- lScore + (sqrt(varLscore)*qnorm(quantile))

  # contrast estimates and confidence bounds
  ve <- VE(markGrid, alphaHat, betaHat, gammaHat)
  if (contrast=="ve"){
    est <- ve
    se <- seVE(markGrid, object$cov, alphaHat, betaHat, gammaHat)
    lb <- c(est - qnorm(quantile)*se)
    ub <- c(est + qnorm(quantile)*se)
  } else if (contrast=="hr"){
    est <- 1 - ve
    lb <- exp(lLB)
    ub <- exp(lUB)
  } else if (contrast=="loghr"){
    est <- log(1 - ve)
    lb <- lLB
    ub <- lUB
  }


  # create contrast matrix to be included in output
  contrastMatrix <- cbind(markGrid, est, lb, ub)
  contrastName <- ifelse(contrast=="ve", "VE", ifelse(contrast=="hr", "HR", "LogHR"))
  colnames(contrastMatrix) <- c(colnames(object$mark), contrastName, "LB", "UB")


  out <- list(coef, pLR.HRunity.2sided, pWald.HRunity.2sided, pWtWald.HRunity.1sided,
               pWald.HRconstant, pLR.HRconstant, contrastMatrix)
  names(out) <- c("coef", "pLR.HRunity.2sided", "pWald.HRunity.2sided", "pWtWald.HRunity.1sided",
                  paste0("pWald.HRconstant", ifelse(sieveAlternative=="twoSided", "2sided", "1sided")),
                  paste0("pLR.HRconstant", ifelse(sieveAlternative=="twoSided", "2sided", "1sided")),
                  "contrast")

  class(out) <- "summary.sievePH"
  return(out)
}
