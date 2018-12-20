dVE <- function(v, alpha, beta, gamma){ -exp(alpha+beta*v+gamma)*c(1,v,1) }

# 'seVE' calculates the standard error given parameters for the grid of mark values, v,
# values for alpha, beta, and gamma, and the corresponding covariance matrix
seVE <- function(v, cov, alpha, beta, gamma){
  sapply(v, function(mark){ drop(sqrt(t(dVE(mark,alpha,beta,gamma)) %*% cov %*% dVE(mark,alpha,beta,gamma))) })
}

LRtest <- function(mark, txInd, thetaHat, lambdaHat){
  V <- cbind(1,mark)
  z <- txInd
  nmark <- NCOL(V)
  
  g <- function(theta){ exp(drop(V %*% theta)) }
  loglik <- function(theta, lambda){ -sum(log(1 + lambda*(g(theta)-1))) + sum(z*log(g(theta))) }
  
  teststat <- 2*(loglik(thetaHat, lambdaHat) - loglik(rep(0,nmark), 0))
  pval <- 1-pchisq(teststat,nmark-1)
  list(teststat=teststat, pval=pval)
}

#' Summarizing Results for Semiparametric Efficient Estimation and Testing for a Multivariate Mark-Specific Hazard Ratio Model
#' 
#' summary method for an object of class "sievePH"
#' 
#' @param object an object of class "sievePH", usually a result of a call to \code{\link{sievePH}}.
#' @param mark a numeric matrix specifying the values of the multivariate mark variable, where rows correspond to subjects and columns correspond to components of the mark.
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
#' \item \code{coef}{a numeric matrix.}
#' \item \code{pLR.HRunity.2sided}{the p-value for the two-sided likelihood ratio test of unity.}
#' \item \code{pWald.HRunity.2sided}{the p-value for the two-sided WAld test of unity.}
#' \item \code{pWtWald.HRunity.1sided}{the p-value for the one-sided weighted Wald test of unity.}
#' \item \code{pLR.HRconstant.2sided}{the p-value for the two-sided likelihood ratio test of constancy.}
#' \item \code{pWald.HRconstant.2sided}{the p-value for the two-sided Wald test of constancy.}
#' \item \code{ve}{a numeric matrix containing the values of the multivariate mark variable and point and interval estimates for the component specified by the \code{contrast} input parameter, with the confidence level specified by the \code{confLevel} input parameter.}
#'    
#' @examples
#' ### from the example of 'sievePH':   
#' fit <- sievePH(...)
#' 
#' summary(fit)     

summary.sievePH <- function(object, mark, 
                            contrast = c("ve", "hr", "loghr"), 
                            sieveAlternative = c("twoSided","oneSided"), confLevel = 0.95) {
  
  sieveAlternative <- match.arg(sieveAlternative, choices = c("twoSided","oneSided"))
  contrast <- match.arg(contrast, choices = c("ve", "hr", "loghr"))
  
  betaHat <- object$betaHat
  alphaHat <- object$alphaHat
  lambdaHat <- object$lambdaHat
  gammaHat <- object$gammaHat
  vAlphaHat <- object$cov[1,1]
  vBetaHat <- object$cov[2,2]
  vGammaHat <- object$cov[3,3]
  V <- object$mark
  Z <- object$txInd
  ve <- object$ve
  
  # standard error and confidence interval limits
  se <- seVE(V,cov,alphaHat,betaHat,gammaHat)
  quantile <- 1 - (1-confLevel)/2
  lb <- c(ve - qnorm(quantile)*se)
  ub <- c(ve + qnorm(quantile)*se)
  
  ### two-sided likelihood ratio test of unity 
  ### H00: beta=0 and gamma=0 (i.e., H00: VE(v)=0) vs. H1: beta!=0 or gamma!=0
  
  ### two-sided Wald test of unity 
  ### H00: beta=0 and gamma=0 (i.e., H00: VE(v)=0) vs. H1: beta!=0 or gamma!=0
  waldH00 <- drop(t(c(thetaHat[2], gammaHat)) %*% solve(Sigma[2:3,2:3]) %*% c(thetaHat[2], gammaHat))
  pWald.HR.unity.2sided <- 2*(1 - pnorm(abs(waldH00)))
  
  ### one-sided weighted Wald-type test of unity 
  ### H00: HR(v)=1 vs H1: HR<1 are HR(v) increasing in each component of v
  weighted.waldH00 <- (betaHat/vBetaHat - gammaHat/vGammaHat)/
    sqrt(1/vBetaHat + 1/vGammaHat - 2*cov[3,2]/(vBetaHat*vGammaHat))
  pWtWald.HRunity.1sided <- 1 - pnorm(weighted.waldH00)
  
  ### two-sided Wald test of H0: VE(v)=VE
  waldH0beta <- betaHat/sqrt(vBetaHat)
  waldH0beta.pval <- 2*(1 - pnorm(abs(waldHObeta)))
  
  ### two-sided Wald test of H0: alpha=0
  waldH0alpha <- alphaHat/sqrt(vAlphaHat)
  waldH0alpha.pval <- 2*(1 - pnorm(abs(waldH0alpha)))
  
  ### two-sided Wald test of H0: gamma=0
  waldH0gamma <- gammaHat/sqrt(vGammaHat)
  waldH0gamma.pval <- 2*(1 - pnorm(abs(waldH0gamma)))
  
  if (sieveAlternative == "oneSided") {
    
    ### 1-sided Wald test of H0:VE(v)=VE (beta=0) vs alternative that beta > 0
    waldH0 <- betaHat/sqrt(vBetaHat)
    pWald.HRconstant <- 1 - pnorm(waldH0)
    
    ### 1-sided likelihood ratio test of H0: VE(v)=VE (beta=0) vs alternative that beta > 0 
    pLR.HRconstant <- LRtest(object$mark[d==1], object$txInd[d==1], c(alphaHat, betaHat), lambdaHat)$pval
    
  } else {
    
    ### 2-sided Wald test of H0:VE(v)=VE (beta=0)
    waldH0 <- betaHat/sqrt(vBetaHat)
    pWald.HRconstant <- 2*(1 - pnorm(abs(waldH0)))
    
    ### 2-sided likelihood ratio test of H0: VE(v)=VE (beta=0) 
    pLR.HRconstant <- LRtest(object$mark[d==1], object$txInd[d==1], c(alphaHat, betaHat), lambdaHat)$pval
    
  }
  
  # tab <- cbind(ve, se, lb, ub, waldH0beta.pval, waldH0alpha.pval, waldH0gamma.pval, 
  #              weighted.waldH00.pval, waldH0.pval, lrBeta.pval)
  # colnames(tab) <- c("VE","SE","LB","UB","pWaldH0beta", "pWaldH0alpha", "pWaldH0gamma", 
  #                    "pWeightedWaldH00", paste0("pWaldH0", ifelse(oneSided, "1sided", "2sided")), 
  #                    paste0("pLRH0", ifelse(oneSided, "1sided", "2sided")))
  
  coef <- matrix()
  
  ve <- cbind(mark, ve, lb, ub)
  contrastName <- ifelse(contrast=="ve", "VE", ifelse(contrast=="hr", "HR", "LogHR"))
  colnames(ve) <- c(colnames(mark), contrastName, "LB", "UB")
  
  
  out <- list(coef, pLR.HRunity.2sided, pWald.HRunity.2sided, pWtWald.HRunity.1sided,
               pWald.HRconstant, pLR.HRconstant, ve)
  names(out) <- c("coef", "pLR.HRunity.2sided", "pWald.HRunity.2sided", "pWtWald.HRunity.1sided", 
                  paste0("pWald.HRconstant", ifelse(sieveAlternative=="twoSided", "2sided", "1sided")), 
                  paste0("pLR.HRconstant", ifelse(sieveAlternative=="twoSided", "2sided", "1sided")),
                  contrast)
  
  class(out) <- "summary.sievePH"
  return(out)
}
