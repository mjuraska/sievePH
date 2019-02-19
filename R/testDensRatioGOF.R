### goodness-of-fit test of the validity of a density ratio model (Qin and Zhang, 1997)
### can handle univariate and bivariate marks
### B = number of bootstrap iterations

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

#' Goodness-of-Fit Test of the Validity of a Univariate or Multivariate Mark Density Ratio Model
#'
#' \code{testDensRatioGoF} implements the complete-case goodness-of-fit test of Qin and Zhang (1997) for the validity of the mark density ratio model in 
#' Juraska and Gilbert (2013). \code{testDensRatioGoF} can handle univariate and multivariate marks, but the mark data must be fully observed in all failures.
#'
#' @param mark either a numeric vector specifying a univariate continuous mark or a matrix specifying a multivariate continuous mark, 
#' completely observed in all cases (i.e., failures). No missing mark values are permitted.
#' @param tx a numeric vector indicating the treatment group (1 if treatment, 0 if placebo) in all cases. No missing values are permitted.
#' @param theta a numeric vector of the coefficients \eqn{\phi} in the weight function \eqn{g(v, \phi)} in the density ratio model
#' @param lambda the Lagrange multiplier in the profile score functions for \eqn{\phi} (that arises by profiling out the nuisance parameter)
#' @param B the number of bootstrap iterations (1000 by default)
#'
#' @return Returns a list containing the following components:
#' \itemize{
#' \item \code{teststat}: the Kolmogorov-Smirnov-type test statistic from the test of validity of the mark density ratio model
#' \item \code{pval}: the bootstrap p-value from the test of validity of the mark density ratio model
#' \item \code{theta}: the input object, or a numeric vector of estimates of coefficients \eqn{\phi} in the weight function \eqn{g(v, \phi)} in the density ratio model
#' \item \code{lambda}: the input object, or an estimate of the Lagrange multiplier in the profile score functions for \eqn{\phi} (that arises by profiling out the nuisance parameter)
#'
#' @references Juraska, M. and Gilbert, P. B. (2013), Mark-specific hazard ratio model with multivariate continuous marks: an application to vaccine efficacy. \emph{Biometrics} 69(2):328-337.
#'
#' Qin, J., & Zhang, B. (1997). A goodness-of-fit test for logistic regression models based on case-control data. \emph{Biometrika}, 84(3), 609-618.
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
#' mark1 <- mark1[eventInd==1]
#' mark2 <- mark2[eventInd==1]
#' tx <- tx[eventInd==1]
#'
#' # perform the test with a univariate mark
#' testDensRatioGOF(mark1, tx, B=20)
#' 
#' # perform the test with a bivariate mark
#' testDensRatioGOF(cbind(mark1, mark2), tx, B=20)
#'
#' @export
testDensRatioGOF <- function(mark, tx, theta=NULL, lambda=NULL, B=1000){
  ninf <- length(tx)   ## number of infections
  n0 <- sum(1-tx)      ## number of placebo infections
  n1 <- sum(tx)        ## number of vaccine infections
  if (sum(is.null(theta),is.null(lambda))>0){
    param <- densRatio(mark, tx)$coef
    theta <- param[-length(param)]
    lambda <- param[length(param)]
  }
  
  g <- function(mark, theta){ exp(drop(cbind(1,mark) %*% theta)) }
  p <- function(mark, theta, lambda, m){ 1/(m*(1+lambda*(g(mark,theta)-1))) }
  F0np <- function(mark, tx){
    if (NCOL(as.matrix(mark))==1){
      F0np.vector <- sapply(mark, function(mark.i){ mean(mark[tx==0]<=mark.i) })
    } else {
      F0np.vector <- apply(mark, 1, function(mark.i){
        mark.vector <- sapply(1:NCOL(as.matrix(mark)), function(col){ mark[tx==0, col] <= mark.i[col] })
        mean(Reduce("&", as.list(as.data.frame(mark.vector)))) 
      })
    }
    return( F0np.vector )
  }
  F0sp <- function(mark, prob){
    if (NCOL(as.matrix(mark))==1){
      F0sp.vector <- sapply(mark, function(mark.i){ sum(prob*ifelse(mark<=mark.i,1,0)) })
    } else { 
      F0sp.vector <- apply(mark, 1, function(mark.i){ 
        mark.vector <- sapply(1:NCOL(as.matrix(mark)), function(col){ mark[, col] <= mark.i[col] })
        sum(prob*ifelse(Reduce("&", as.list(as.data.frame(mark.vector))),1,0)) 
      })
    }
    return( F0sp.vector )
  }
  
  f0 <- p(mark,theta,lambda,ninf)
  delta <- sqrt(ninf)*max(abs(F0np(mark,tx) - F0sp(mark,f0)))
  
  N0.ast <- rmultinom(B,n0,f0)
  N1.ast <- rmultinom(B,n1,f0*g(mark,theta))
  v0.ast <- lapply(1:B, function(iter){
    if (NCOL(as.matrix(mark))==1){
      out <- rep(mark,N0.ast[,iter])
    } else {
      out <- cbind(sapply(1:NCOL(as.matrix(mark)), function(col){ rep(mark[, col], N0.ast[, iter])}))
    }
    return( out )
  })
  v1.ast <- lapply(1:B, function(iter){
    if (NCOL(as.matrix(mark))==1){
      out <- rep(mark,N1.ast[,iter])
    } else {
      out <- cbind(sapply(1:NCOL(as.matrix(mark)), function(col){ rep(mark[, col], N1.ast[, iter])}))
    }
    return( out )
  })
  z.ast <- rep(0:1,c(n0,n1))
  
  teststat <- sapply(1:B, function(iter){
    v.ast.iter <- rbind(as.matrix(v0.ast[[iter]]),as.matrix(v1.ast[[iter]]))
    param <- densRatio(v.ast.iter, z.ast)$coef
    theta.ast <- param[-length(param)]
    lambda.ast <- param[length(param)]
    sqrt(ninf)*max(abs(F0np(v.ast.iter,z.ast) - F0sp(v.ast.iter,p(v.ast.iter,theta.ast,lambda.ast,ninf))))
  })
  pval <- mean(teststat>=delta)
  return(list(teststat=delta, pval=pval, theta=theta, lambda=lambda))
}
