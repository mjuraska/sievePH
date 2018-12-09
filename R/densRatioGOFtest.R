### goodness-of-fit test of the validity of a density ratio model (Qin and Zhang, 1997)
### can handle univariate and bivariate marks
### B = number of bootstrap iterations
densRatioGOFtest <- function(mark, trt.id, theta=NULL, lambda=NULL, B=1000){
  ninf <- length(trt.id)   ## number of infections
  n0 <- sum(1-trt.id)      ## number of placebo infections
  n1 <- sum(trt.id)        ## number of vaccine infections
  if (sum(is.null(theta),is.null(lambda))>0){
    param <- densRatio(mark, trt.id)$coef
    theta <- param[-length(param)]
    lambda <- param[length(param)]
  }
  
  g <- function(mark, theta){ exp(drop(cbind(1,mark) %*% theta)) }
  p <- function(mark, theta, lambda, m){ 1/(m*(1+lambda*(g(mark,theta)-1))) }
  F0np <- function(mark, trt.id){
    if (NCOL(as.matrix(mark))==1){
      F0np.vector <- sapply(mark, function(mark.i){ mean(mark[trt.id==0]<=mark.i) })
    } else {
      F0np.vector <- apply(mark, 1, function(mark.i){ mean(mark[trt.id==0,1]<=mark.i[1] & mark[trt.id==0,2]<=mark.i[2]) })
    }
    return( F0np.vector )
  }
  F0sp <- function(mark, prob){
    if (NCOL(as.matrix(mark))==1){
      F0sp.vector <- sapply(mark, function(mark.i){ sum(prob*ifelse(mark<=mark.i,1,0)) })
    } else {
      F0sp.vector <- apply(mark, 1, function(mark.i){ sum(prob*ifelse(mark[,1]<=mark.i[1] & mark[,2]<=mark.i[2],1,0)) })
    }
    return( F0sp.vector )
  }
  
  f0 <- p(mark,theta,lambda,ninf)
  delta <- sqrt(ninf)*max(abs(F0np(mark,trt.id) - F0sp(mark,f0)))
  
  N0.ast <- rmultinom(B,n0,f0)
  N1.ast <- rmultinom(B,n1,f0*g(mark,theta))
  v0.ast <- lapply(1:B, function(iter){
    if (NCOL(as.matrix(mark))==1){
      out <- rep(mark,N0.ast[,iter])
    } else {
      out <- cbind(rep(mark[,1],N0.ast[,iter]), rep(mark[,2],N0.ast[,iter]))
    }
    return( out )
  })
  v1.ast <- lapply(1:B, function(iter){
    if (NCOL(as.matrix(mark))==1){
      out <- rep(mark,N1.ast[,iter])
    } else {
      out <- cbind(rep(mark[,1],N1.ast[,iter]), rep(mark[,2],N1.ast[,iter]))
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
  list(teststat=delta, pval=pval, theta=theta, lambda=lambda)
}