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

densRatioAUG <- function(mark, trt.id, aux.miss, aux){
	V <- cbind(1,mark)
	V.complete <- na.omit(V)
	z <- trt.id
	na.idx <- attr(V.complete,"na.action")
	if (!is.null(na.idx)){
		z.complete <- z[-na.idx]
		aux.complete <- aux[-na.idx]
	} else {
		z.complete <- z
		aux.complete <- aux
	}
	nmark <- NCOL(V.complete)
	ninf <- NROW(V.complete)
	VV.complete <- apply(V.complete,1,tcrossprod)
		
	g <- function(theta){ exp(drop(V.complete %*% theta)) }
	dG <- function(theta){ t(g(theta) * V.complete) }
	d2G <- function(theta){ array(t(t(VV.complete)*g(theta)),dim=c(nmark,nmark,ninf)) }
	dGdG <- function(theta){ array(apply(dG(theta),2,tcrossprod),dim=c(nmark,nmark,ninf)) }

	fscore.i1 <- function(theta, lambda){
		-lambda * t(dG(theta)) * (1/(1+lambda*(g(theta)-1))) + t(dG(theta)) * (z.complete/g(theta))
	}
	fscore.i2 <- function(theta, lambda){
		-(g(theta)-1)/(1+lambda*(g(theta)-1))
	}
	fscore.i <- function(theta, lambda){
		cbind(fscore.i1(theta, lambda), fscore.i2(theta, lambda))
	}
	aug.mean1 <- function(theta, lambda){
		U <- fscore.i1(theta, lambda)
		predicted.vals <- sapply(1:NCOL(U), function(col){
			fit <- lm(U[,col] ~ z.complete*aux.complete + I(aux.complete^2))
			predict(fit, data.frame(z.complete=z, aux.complete=aux))
		})
		return( predicted.vals )
	}
	aug.mean2 <- function(theta, lambda){
		U <- fscore.i2(theta, lambda)
		fit <- lm(U ~ z.complete*aux.complete + I(aux.complete^2))
		predict(fit, data.frame(z.complete=z, aux.complete=aux))
	}

	score1 <- function(theta, lambda){
		drop(-lambda * dG(theta) %*% (1/(pi*(1+lambda*(g(theta)-1)))) + dG(theta) %*% (z.complete/(pi*g(theta))) +
		t(aug.mean1(theta, lambda)) %*% (1-r/pi.all))
	}
	score2 <- function(theta, lambda){
		-sum((g(theta)-1)/(pi*(1+lambda*(g(theta)-1)))) + sum(aug.mean2(theta, lambda)*(1-r/pi.all))
	}
	score <- function(theta, lambda){ c(score1(theta,lambda), score2(theta,lambda)) }
	jack11 <- function(theta, lambda){
		d2Gperm <- aperm(d2G(theta), c(3,1,2))
		dGdGperm <- aperm(dGdG(theta), c(3,1,2))
		term1 <- apply(aperm(d2Gperm*(1/(pi*(1+lambda*(g(theta)-1)))), c(2,3,1)),c(1,2),sum)
		term2 <- apply(aperm(dGdGperm*(1/(pi*(1+lambda*(g(theta)-1))^2)), c(2,3,1)),c(1,2),sum)
		term3 <- apply(aperm(d2Gperm*(z.complete/(pi*g(theta))), c(2,3,1)),c(1,2),sum)
		term4 <- apply(aperm(dGdGperm*(z.complete/(pi*g(theta)^2)), c(2,3,1)),c(1,2),sum)

		d2U.theta1 <- aperm(d2Gperm*(1/(1+lambda*(g(theta)-1))), c(2,3,1))
		d2U.theta2 <- aperm(dGdGperm*(1/(1+lambda*(g(theta)-1))^2), c(2,3,1))
		d2U.theta3 <- aperm(d2Gperm*(z.complete/g(theta)), c(2,3,1))
		d2U.theta4 <- aperm(dGdGperm*(z.complete/g(theta)^2), c(2,3,1))
		d2U.theta <- -lambda*(d2U.theta1 - lambda*d2U.theta2) + d2U.theta3 - d2U.theta4
		predicted.vals.jack <- array(0, dim=c(nmark,nmark,length(z)))
		for (i in 1:nmark){
		for (j in 1:nmark){
			resp <- d2U.theta[i,j,]
			fit <- lm(resp ~ z.complete*aux.complete + I(aux.complete^2))
			predicted.vals.jack[i,j,] <- predict(fit, data.frame(z.complete=z, aux.complete=aux))
		}}
		weighted.predicted.vals.jack <- apply(aperm(aperm(predicted.vals.jack, c(3,1,2))*(1-r/pi.all), c(2,3,1)), c(1,2), sum)

		return( -lambda*(term1 - lambda*term2) + term3 - term4 + weighted.predicted.vals.jack )
	}
	jack21 <- function(theta, lambda){
		d2U.theta.lambda <- -t(dG(theta)) * (1/(1+lambda*(g(theta)-1))^2)
		predicted.vals.jack <- matrix(0,nrow=length(z),ncol=nmark)
		for (i in 1:nmark){
			resp <- d2U.theta.lambda[,i]
			fit <- lm(resp ~ z.complete*aux.complete + I(aux.complete^2))
			predicted.vals.jack[,i] <- predict(fit, data.frame(z.complete=z, aux.complete=aux))
		}
		weighted.predicted.vals.jack <- colSums(predicted.vals.jack*(1-r/pi.all))
		return( drop(-dG(theta) %*% (1/(pi*(1+lambda*(g(theta)-1))^2))) + weighted.predicted.vals.jack )
	}
	jack22 <- function(theta, lambda){
		d2U.lambda <- ((g(theta)-1)/(1+lambda*(g(theta)-1)))^2
		fit <- lm(d2U.lambda ~ z.complete*aux.complete + I(aux.complete^2))
		predicted.vals.jack <- predict(fit, data.frame(z.complete=z, aux.complete=aux))
		weighted.predicted.vals.jack <- sum(predicted.vals.jack*(1-r/pi.all))
		return( sum(((g(theta)-1)^2)/(pi*(1+lambda*(g(theta)-1))^2)) + weighted.predicted.vals.jack )
	}
	jack <- function(theta, lambda){
		j21 <- jack21(theta,lambda)
		cbind(rbind(jack11(theta,lambda),j21),c(j21,jack22(theta,lambda)))
	}

	r <- apply(V, 1, function(row){ ifelse(sum(is.na(row))>0,0,1) })
	# pi.all <- glm(r ~ z + aux.miss + z*aux.miss, family=binomial)$fitted
	pi.all <- glm(r ~ z, family=binomial)$fitted
	if (!is.null(na.idx)){
		pi <- pi.all[-na.idx]
	} else {
		pi <- pi.all
	}
	if (any(pi==0)){ stop("Selection probabilities not bounded away from 0.") }

	param.old <- numeric(nmark+1)
	param.new <- c(numeric(nmark),0.5)
#	param.old <- densRatioIPW(mark, trt.id, aux)$coef
#	param.new <- param.old + c(numeric(nmark),1e-2)
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
		U[,1:nmark] <- U[,1:nmark] + aug.mean1(theta, lambda) * (1-r/pi.all)
		U[-na.idx,nmark+1] <- (g(theta)-1)/(pi*(1+lambda*(g(theta)-1)))
		U[,nmark+1] <- U[,nmark+1] + aug.mean2(theta, lambda) * (1-r/pi.all)
		# S <- (r-pi.all) * cbind(1,z,aux.miss,z*aux.miss)
		S <- (r-pi.all) * cbind(1,z)
		# resids <- sapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2] + S[,3] + S[,4])$resid })
		resids <- sapply(1:NCOL(U), function(i){ lm(U[,i] ~ S[,1] + S[,2])$resid })
		crossprod(resids)/ninf
	}

	JackInv <- try(solve(jack(theta.new,lambda.new)), silent=TRUE)
	if (class(JackInv)!="try-error"){
		Var <- ninf * JackInv %*% Resid(theta.new,lambda.new) %*% JackInv
		names(param.new) <- rownames(Var) <- colnames(Var) <- c("alpha",paste("beta",1:(nmark-1),sep=""),"lambda")
	} else {
		Var <- NULL
	}

	list(coef=param.new, var=Var, jack=jack11(theta.new,lambda.new), conv=!(class(jackInv)=="try-error" | class(JackInv)=="try-error"))
}

####################################################################################
### Description: R functions to perform the test of independence between T and V ###
###              based on the Huang-Louis estimator for the joint cdf of (T,V)   ###
###              (Huang and Louis, Biometrika, 1998)                             ###
###              (the most recent update of 'ftu')                               ###
### Author:      Michal Juraska                                                  ###
### Date:        November 19, 2011                                               ###
####################################################################################

# T: survival time
# U: mark
# C: censoring time
# data: replicates of [X=min(T,C), Delta=I(T<=C), Y=U*Delta]

# input: data[size,3] each row corresponds to a data point (X,Delta,Y)
#        tu[,2] each row corresponds to the coordinates (t,u) of a point
# output: a vector of the same length as tu, cdf values at (t,u)

ftu <- function(tu,data){
  ord <- order(data[,1],-data[,2])
  size <- length(data[,1])
  prob <- numeric(size)
  surv <- 1
  for (i in 1:size){
    if (data[ord[i],2]==1){ prob[ord[i]] <- surv/(size-i+1) }
    surv <- surv-prob[ord[i]]
  }
  sapply(1:NROW(tu), function(j){
    sum(prob[data[,1]<=tu[j,1] & data[,3]<=tu[j,2]],na.rm=TRUE)
  })
}

g <- function(tu,data){
  if (is.vector(tu)) tu <- t(as.matrix(tu)) 
  tu <- na.omit(tu)
  abs(ftu(tu,data) - ftu(cbind(tu[,1],Inf),data)*ftu(cbind(Inf,tu[,2]),data[data[,2]==1,]))
}

### bootstrap to estimate the p-value for the test of independence between T and V
bootpval <- function(data, iter=1000){
  T <- max(g(data[data[,2]==1,c(1,3)],data))
  n <- NROW(data)
  resamp <- matrix(sample(1:n,n*iter,replace=TRUE),n,iter)
  obs.mark <- as.vector(data[data[,2]==1,3])
  bT <- sapply(1:iter, function(j){
    bdata <- data[resamp[,j],1:2]
    if (sum(bdata[,2])>1){
      m <- sum(bdata[,2])
      bmark <- numeric(n)
      bmark[which(bdata[,2]==1)] <- sample(obs.mark,m,replace=TRUE)
      bdata <- cbind(bdata,bmark)
      tstat <- max(g(bdata[bdata[,2]==1,c(1,3)],bdata))
    } else {
      tstat <- NULL
    }
    tstat
  })
  if (is.list(bT)){ bT <- do.call("c",lapply(bT,"[",1)) }
  mean(bT>=T)
}

### estimator for \cov(\hat{\phi},\hat{\lambda}) using Theorem 1 in Juraska and Gilbert (2013, Biometrics)

library(survival)

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
  omega <- drop(score1.vect(phi.hat,lambda.hat) %*% (score3.vect(gamma.hat) + p*l.vect(gamma.hat))/n - 
                  sum(score3.vect(gamma.hat) + p*l.vect(gamma.hat))*apply(score1.vect(phi.hat,lambda.hat),1,sum)/(n^2))   # a vector with 2 components
  drop(solve(jack(phi.hat,lambda.hat))[1:2,1:2] %*% omega)/(n*jack33)
}

### estimator for \cov(\hat{\phi}_{ipw},\hat{\lambda}) based on Theorem 1 of JG (2013)
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

### estimator for \cov(\hat{\phi}_{aug},\hat{\lambda}) based on Theorem 1 of JG (2013)
covEstAUG <- function(time, find, mark, txind, aux.miss, aux, phi.hat, lambda.hat, gamma.hat){
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
  aux.complete <- aux[find==1 & !is.na(mark)]
  aux.f <- aux[find==1]
  aux.miss.f <- aux.miss[find==1]
  time.completeM <- matrix(time.complete, nrow=n, ncol=m.complete, byrow=TRUE)
  VV.complete <- apply(V.complete,1,tcrossprod)
  nmark <- NCOL(V.complete)
  
  g <- function(phi){ exp(drop(V.complete %*% phi)) }
  dG <- function(phi){ t(g(phi) * V.complete) }
  d2G <- function(phi){ array(t(t(VV.complete)*g(phi)),dim=c(nmark,nmark,m.complete)) }
  dGdG <- function(phi){ array(apply(dG(phi),2,tcrossprod),dim=c(nmark,nmark,m.complete)) }
  
  score1.complete.vect <- function(phi, lambda){
    (-lambda/(1+lambda*(g(phi)-1)) + txind.complete/g(phi)) * t(dG(phi))
  }
  score2.complete.vect <- function(phi, lambda){
    -(g(phi)-1)/(1+lambda*(g(phi)-1))
  }
  aug.mean1 <- function(phi, lambda){
    U <- score1.complete.vect(phi, lambda)
    predicted.vals <- sapply(1:NCOL(U), function(col){
      fit <- lm(U[,col] ~ txind.complete*aux.complete + I(aux.complete^2))
      predict(fit, data.frame(txind.complete=txind.f, aux.complete=aux.f))
    })
    predicted.vals
  }
  aug.mean2 <- function(phi, lambda){
    U <- score2.complete.vect(phi, lambda)
    fit <- lm(U ~ txind.complete*aux.complete + I(aux.complete^2))
    predict(fit, data.frame(txind.complete=txind.f, aux.complete=aux.f))
  }
  score1.vect <- function(phi, lambda){
    vect <- matrix(0, nrow=m.f, ncol=nmark)
    vect[-na.idx,] <- (-lambda/(pi*(1+lambda*(g(phi)-1))) + txind.complete/(pi*g(phi))) * t(dG(phi))
    vect <- vect + (aug.mean1(phi, lambda) * (1-r/pi.all))
    t(vect)
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
    drop(-lambda * dG(phi) %*% (1/(pi*(1+lambda*(g(phi)-1)))) + dG(phi) %*% (txind.complete/(pi*g(phi))) +
           t(aug.mean1(phi, lambda)) %*% (1-r/pi.all))
  }
  score2 <- function(phi, lambda){
    -sum((g(phi)-1)/(pi*(1+lambda*(g(phi)-1)))) + sum(aug.mean2(phi, lambda)*(1-r/pi.all))
  }
  score <- function(phi, lambda){ c(score1(phi,lambda),score2(phi,lambda)) }
  jack11 <- function(phi, lambda){
    d2Gperm <- aperm(d2G(phi), c(3,1,2))
    dGdGperm <- aperm(dGdG(phi), c(3,1,2))
    term1 <- apply(aperm(d2Gperm*(1/(pi*(1+lambda*(g(phi)-1)))), c(2,3,1)),c(1,2),sum)
    term2 <- apply(aperm(dGdGperm*(1/(pi*(1+lambda*(g(phi)-1))^2)), c(2,3,1)),c(1,2),sum)
    term3 <- apply(aperm(d2Gperm*(txind.complete/(pi*g(phi))), c(2,3,1)),c(1,2),sum)
    term4 <- apply(aperm(dGdGperm*(txind.complete/(pi*g(phi)^2)), c(2,3,1)),c(1,2),sum)
    
    d2U.phi1 <- aperm(d2Gperm*(1/(1+lambda*(g(phi)-1))), c(2,3,1))
    d2U.phi2 <- aperm(dGdGperm*(1/(1+lambda*(g(phi)-1))^2), c(2,3,1))
    d2U.phi3 <- aperm(d2Gperm*(txind.complete/g(phi)), c(2,3,1))
    d2U.phi4 <- aperm(dGdGperm*(txind.complete/g(phi)^2), c(2,3,1))
    d2U.phi <- -lambda*(d2U.phi1 - lambda*d2U.phi2) + d2U.phi3 - d2U.phi4
    predicted.vals.jack <- array(0, dim=c(nmark,nmark,m.f))
    for (i in 1:nmark){
      for (j in 1:nmark){
        resp <- d2U.phi[i,j,]
        fit <- lm(resp ~ txind.complete*aux.complete + I(aux.complete^2))
        predicted.vals.jack[i,j,] <- predict(fit, data.frame(txind.complete=txind.f, aux.complete=aux.f))
      }}
    weighted.predicted.vals.jack <- apply(aperm(aperm(predicted.vals.jack, c(3,1,2))*(1-r/pi.all), c(2,3,1)), c(1,2), sum)
    -lambda*(term1 - lambda*term2) + term3 - term4 + weighted.predicted.vals.jack
  }
  jack21 <- function(phi, lambda){
    d2U.phi.lambda <- -t(dG(phi)) * (1/(1+lambda*(g(phi)-1))^2)
    predicted.vals.jack <- matrix(0,nrow=m.f,ncol=nmark)
    for (i in 1:nmark){
      resp <- d2U.phi.lambda[,i]
      fit <- lm(resp ~ txind.complete*aux.complete + I(aux.complete^2))
      predicted.vals.jack[,i] <- predict(fit, data.frame(txind.complete=txind.f, aux.complete=aux.f))
    }
    weighted.predicted.vals.jack <- colSums(predicted.vals.jack*(1-r/pi.all))
    drop(-dG(phi) %*% (1/(pi*(1+lambda*(g(phi)-1))^2))) + weighted.predicted.vals.jack
  }
  jack22 <- function(phi, lambda){
    d2U.lambda <- ((g(phi)-1)/(1+lambda*(g(phi)-1)))^2
    fit <- lm(d2U.lambda ~ txind.complete*aux.complete + I(aux.complete^2))
    predicted.vals.jack <- predict(fit, data.frame(txind.complete=txind.f, aux.complete=aux.f))
    weighted.predicted.vals.jack <- sum(predicted.vals.jack*(1-r/pi.all))
    sum(((g(phi)-1)^2)/(pi*(1+lambda*(g(phi)-1))^2)) + weighted.predicted.vals.jack
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

### Usage:
# covEst(X,d,V,Z,dRatio$coef[1:2],dRatio$coef[3],gammaHat)
# covEstIPW(X,d,V,Z,A.miss,dRatio$coef[1:2],dRatio$coef[3],gammaHat)
# covEstAUG(X,d,V,Z,A.miss,A.miss,dRatio$coef[1:2],dRatio$coef[3],gammaHat)

