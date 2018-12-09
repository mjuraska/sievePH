seVE <- function(v, Var, a, b, g){
  sapply(v, function(mark){ drop(sqrt(t(dVE(mark,a,b,g)) %*% Var %*% dVE(mark,a,b,g))) })
}

LRtest <- function(mark, trt.id, theta.hat, lambda.hat){
  V <- cbind(1,mark)
  z <- trt.id
  nmark <- NCOL(V)
  
  g <- function(theta){ exp(drop(V %*% theta)) }
  loglik <- function(theta, lambda){ -sum(log(1 + lambda*(g(theta)-1))) + sum(z*log(g(theta))) }
  
  teststat <- 2*(loglik(theta.hat, lambda.hat) - loglik(rep(0,nmark), 0))
  pval <- 1-pchisq(teststat,nmark-1)
  list(teststat=teststat, pval=pval)
}

#' Summarizing Results for Semiparametric Efficient Estimation and Testing for a Multivariate Mark-Specific Hazard Ratio Model
#' 
#' summary method for an object of class "sievePH"
#' 
#' @param object an object of class "sievePH", usually a result of a call to \code{\link{sievePH}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @details 
#' \code{print.summary.speff} prints a formatted summary of results. An inferential table is 
#' produced with point and interval estimates for the mark-specific hazard ratio, standard 
#' error estimates, and Wald test p-values for the two relevant hypothesis tests of no 
#' vaccine efficacy and no sieve effect.
#' 
#' @return 
#' A list with the following components:
#' \itemize{
#'   \item \code{tab} an inferential table for the mark-specific hazard ratio
#' }
#'    
#' @examples
#' ### from the example of 'sievePH':   
#' fit <- sievePH(...)
#' 
#' summary(fit)     

summary.sievePH <- function(object,...) {
  
  betaHat <- object$betaHat
  alphaHat <- object$alphaHat
  lambdaHat <- object$lambdaHat
  gammaHat <- object$gammaHat
  vAlphaHat <- object$Sigma[1,1]
  vBetaHat <- object$Sigma[2,2]
  vGammaHat <- object$Sigma[3,3]
  V <- object$mark
  Z <- object$txInd
  ve <- object$ve
  
  se <- seVE(V,Sigma,alphaHat,betaHat,gammaHat)
  ci <- c(ve + qnorm(0.975)*se %o% c(-1,1))
  
  ### two-sided Wald test of H0: VE(v)=VE
  waldH0beta <- betaHat/sqrt(vBetaHat)
  waldH0beta.pval <- 2*(1 - pnorm(abs(waldHObeta)))
  
  ### two-sided Wald test of H0: alpha=0
  waldH0alpha <- alphaHat/sqrt(vAlphaHat)
  waldH0alpha.pval <- 2*(1 - pnorm(abs(waldH0alpha)))
  
  ### two-sided Wald test of H0: gamma=0
  waldH0gamma <- gammaHat/sqrt(vGammaHat)
  waldH0gamma.pval <- 2*(1 - pnorm(abs(waldH0gamma)))
  
  ### one-sided weighted Wald-type test of H00: VE(v)=0 vs alternatives where VE>0 and VE(v) is decreasing
  weighted.waldH00 <- (betaHat/vBetaHat - gammaHat/vGammaHat)/
    sqrt(1/vBetaHat + 1/vGammaHat - 2*Sigma[3,2]/(vBetaHat*vGammaHat))
  weighted.waldH00.pval <- 1 - pnorm(weighted.waldH00)
  
  if (oneSided) {
    
    # 1-sided test of H0:VE(v)=VE (beta=0) vs alternative that beta > 0
    waldH0 <- betaHat/sqrt(vBetaHat)
    waldH0.pval <- 1 - pnorm(waldH0)
    
    ### 1-sided likelihood ratio test of H0: VE(v)=VE (beta=0) vs alternative that beta > 0 
    lrBeta.pval <- LRtest(object$mark[d==1], object$txInd[d==1], c(alphaHat, betaHat), lambdaHat)$pval
    
  } else {
    
    # 2-sided test of H0:VE(v)=VE (beta=0)
    waldH0 <- betaHat/sqrt(vBetaHat)
    waldH0.pval <- 2*(1 - pnorm(abs(waldH0)))
    
    ### 2-sided likelihood ratio test of H0: VE(v)=VE (beta=0) 
    lrBeta.pval <- LRtest(V[d==1],Z[d==1],thetaHat[-length(thetaHat)],thetaHat[length(thetaHat)])$pval
    
  }

  Tab <- cbind(ve, se, ci[1], ci[2], waldH0beta.pval, waldH0alpha.pval, waldH0gamma.pval, 
               weighted.waldH00.pval, waldH0.pval, lrBeta.pval)
  colnames(Tab) <- c("VE","SE","LB","UB","WaldH0beta p", "WaldH0alpha p", "WaldH0gamma p", 
                     "WeightedWaldH00 p", paste0("WaldH0", ifelse(oneSided, "1sided", "2sided")," p"), 
                     paste0("LRH0", ifelse(oneSided, "1sided", "2sided")," p"))
  out <- list(table=Tab)
  class(out) <- "summary.sievePH"
  out
}
