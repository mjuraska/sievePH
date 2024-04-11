#' @importFrom Rcpp evalCpp
NULL

#' Nonparametric Kernel-Smoothed Stratified Mark-Specific Proportional Hazards Model with a Univariate Continuous Mark, Missing-at-Random in Some Failures
#'
#' \code {cmprskPHContinMark} implements estimation methods of Sun and Gilbert (2012) and hypothesis testing methods of Gilbert and Sun (2015) for a mark-specific
#' proportional hazards model accommodating that some failures have a missing mark. The methods allow separate baseline mark-specific hazard functions
#' for different baseline subgroups. Missing marks are handled via complete-case, inverse-probability weighting (IPW), or augmented IPW (AIPW) approaches.
#'
#' @param eventTime a numeric vector specifying the observed right-censored event time
#' @param eventInd a numeric vector indicating the event of interest (1 if event, 0 if right-censored)
#' @param mark a numeric vector specifying a univariate continuous mark subject to missingness at random. Missing mark values should be set to \code{NA}.
#' For subjects with \code{eventInd = 0}, the value in \code{mark} should also be set to \code{NA}.
#' @param tx a numeric vector indicating the treatment group (1 if treatment, 0 if placebo)
#' @param aux a numeric vector specifying a binary auxiliary covariate predictive of the probability of observing the mark or the value of the mark
#' The mark missingness model only requires that the auxiliary covariates be observed in subjects who experienced the event of interest.
#' For subjects with \code{eventInd = 0}, the value in \code{aux} may be set to \code{NA}. If no auxiliary covariate is used, \code{aux} should be set to the default value of \code{NULL}.
#' @param strata a numeric vector specifying baseline strata (\code{NULL} by default). If specified, a separate mark-specific baseline hazard is assumed for each stratum. If either an IPW or AIPW method is used, a separate logistic regression model for predicting the probability of a missing mark and a separate model for predicting the mark are considered witin each stratum.
#' @param nboot number of bootstrap iterations for simulating the distribution of test statistics for testing \eqn{H_{10}} vs. \eqn{H_{1a}} or \eqn{H_{1m}} and \eqn{H_{20}} vs. \eqn{H_{2a}} or \eqn{H_{2m}}
#' @param missmethod a character string for the estimation procedure to use. Available missing-mark methods include \code{CC}, \code{IPW}, and \code {AIPW}.
#' @param tau a numeric value specifying the end of the follow-up period for conducting the analysis
#' @param tband a numeric value between 0 and \code{tau} specifying the bandwidth of the kernel smoothing function over the failure time
#' @param hband a numeric value between 0 and 1 specifying the bandwidth of the kernel smoothing function over the mark. Larger bandwidths are recommended for higher percentages of missing marks.
#' @param a a numeric value between 0 and 1 specifying the lower bound of the range for testing the null hypothesis \deqn{H_{10}: HR(v) = 1} for v \eqn{\in [a, b]}; \code{a} needs to be \code{1/nvgrid} multiplied by an integer.
#' @param b a numeric value between 0 and 1 specifying the upper bound of the range for testing the null hypothesis \deqn{H_{10}: HR(v) = 1} for v \eqn{\in [a, b]}; \code{b} needs to be \code{1/nvgrid} multiplied by an integer.
#' @param ntgrid an integer value specifying the number of equally spaced time points for which the mark-specific baseline hazard functions are evaluated
#' @param nvgrid an integer value specifying the number of equally spaced mark values between the minimum and maximum of the observed mark values for which the treatment effects are evaluated
#' @param confLevel the confidence level (0.95 by default) of reported confidence intervals
#' @param seed an integer specifying the random number generation seed for reproducing the test statistics and p-values
#'
#' @details
#' \code {cmprskPHContinMark} analyzes data from a randomized placebo-controlled trial that evaluates treatment efficacy for a time-to-event endpoint with a continuous mark.
#' The parameter of interest is the ratio of the conditional mark-specific hazard functions (treatment/placebo), which is based on a stratified mark-specific proportional hazards model.
#' This model assumes no parametric form for the baseline hazard function nor the treatment effect across different mark values. For data with missing marks, two estimation procedures are implemented.
#' The first one uses inverse probability weighting (IPW) of the complete-case estimator, which leverages auxiliary predictors of whether the mark is observed,
#' whereas the second one augments the IPW complete-case estimator with auxiliary predictors of the missing mark value.
#'
#' @return A list containing the following components:
#' \itemize{
#'
#' \item \code{H10}: a data frame with test statistics (first row) and corresponding p-values (second row) for testing \deqn{H_{10}: HR(v) = 1} for v \eqn{\in [a, b]}.
#' The test statistics are calculated for mark values in \eqn{[a1, b]}, where \eqn{a1 = a + (b-a)/nvgrid}.
#' Columns \code{TSUP1} and \code{Tint1} include test statistics and p-values for testing \deqn{H_{10}} vs. \deqn{H_{1a}: HR(v) != 1} for any v \eqn{\in [a, b]} (general alternative).
#' Columns \code{TSUP1m} and \code{Tint1m} include test statistics and p-values for testing \deqn{H_{10}} vs. \deqn{H_{1m}: HR(v) <= 1} with strict inequality for some v \eqn{\in [a, b]} (monotone alternative).
#' Columns \code{TSUP1} and \code{TSUP1m} are for tests based on extensions of Kolmogorov-Smirnov tests.
#' Columns \code{Tint1} and \code{Tint1m} are for tests based on generalizations of Cramer-von Mises tests--integration-based tests.
#'
#' \item \code{H20}: a data frame with test statistics (first row) and corresponding p-values (second row) for testing \deqn{H_{20}}: HR(v) does not depend on v \eqn{\in [a, b]}.
#' The test statistics are calculated for mark values in \eqn{[a2, b]}, where \eqn{a2 = a + 3*(b-a)/nvgrid}.
#' Columns \code{TSUP2} and \code{Tint2} include test statistics and p-values for testing \deqn{H_{20}} vs. \deqn{H_{2a}}: HR depends on v \eqn{\in [a, b]}  (general alternative).
#' Columns \code{TSUP2m} and \code{Tint2m} include test statistics and p-values for testing \deqn{H_{20}} vs. \deqn{H_{2m}}: HR increases as v increases \eqn{\in [a, b]} (monotone alternative).
#' Columns \code{TSUP2} and \code{TSUP2m} are for tests based on extensions of Kolmogorov-Smirnov tests.
#' Columns \code{Tint2} and \code{Tint2m} are for tests based on generalizations of Cramer-von Mises tests--integration-based tests.
#'
#' \item \code{te}: a data frame summarizing point and interval estimates of the mark-specific treatment efficacy on the grid of mark values defined by \code{nvgrid}. The confidence level is specified by \code{confLevel}. \code{te} is based on AIPW estimates.
#' \item \code{estBeta}: a data frame summarizing point estimates and standard errors of the mark-specific coefficients using the complete-case method (\code{betacom} and \code{secom}),
#' the IPW method (\code{betaipw} and \code{seipw}), and the AIPW method (\code{betaaug} and \code{seaug}).
#' \item \code{cBproc1}: a data frame containing mark values in column \code{Mark}, standardized mark values in column \code{Standardized mark},
#' test processes \eqn{Q^{(1)}(v)} for \code{H10} based on observed data sequence in column \code{Observed} and independent sets of normal samples in columns \eqn{S1, S2, \cdots, S_{nboot}}.
#' \item \code{cBproc2}: a data frame containing mark values in column \code{Mark}, standardized mark values in column \code{Standardized mark},
#' test processes \eqn{Q^{(2)}(v)} for \code{H20} based on observed data sequence in column \code{Observed} and independent sets of normal samples in columns \eqn{S1, S2, \cdots, S_{nboot}}.
#' }
#'
#' @references
#' Gilbert, P. B., and Sun, Y. (2015). Inferences on relative failure rates in stratified mark-specific proportional hazards models with missing marks, with application to human immunodeficiency virus vaccine efficacy trials. \emph{Journal of the Royal Statistical Society Series C: Applied Statistics}, 64(1), 49-73.
#'
#' Sun, Y., & Gilbert, P. B. (2012). Estimation of stratified mark‐specific proportional hazards models with missing marks. \emph{Scandinavian Journal of Statistics}, 39(1), 34-52.
#'
#' @examples
#' set.seed(20240410)
#' n <- 200
#' tx <- rep(0:1, each = n / 2)
#' tm <- c(rexp(n / 2, 0.2), rexp(n / 2, 0.2 * exp(-0.4)))
#' cens <- runif(n, 0, 15)
#' eventTime <- pmin(tm, cens, 3)
#' eventInd <- as.numeric(tm <= pmin(cens, 3))
#' mark <- ifelse(eventInd == 1, c(rbeta(n / 2, 2, 1), rbeta(n / 2, 2, 2)), NA)
#' # a binary auxiliary covariate
#' A <- sapply(exp(-0.5 - 0.2 * mark) / (1 + exp(-0.5 - 0.2 * mark)),
#'             function(p){ ifelse(is.na(p), NA, rbinom(1, 1, p)) })
#' linPred <- -0.8 + 0.4 * tx - 0.2 * A
#' probs <- exp(linPred) / (1 + exp(linPred))
#' R <- rep(NA, n)
#' while (sum(R, na.rm = TRUE) < 10){
#'  R[eventInd == 1] <- sapply(probs[eventInd == 1],
#'                             function(p){ rbinom(1, 1, p) })
#' }
#' # a missing-at-random mark
#' mark[eventInd == 1] <- ifelse(R[eventInd == 1] == 1, mark[eventInd == 1], NA)
#' fit <- kernel_sievePH(eventTime, eventInd, mark, tx, A, nboot = 50,
#'                           missmethod = "AIPW", tau = 3, tband = 0.5, hband = 0.3,
#'                           a = 0.1, b = 1, ntgrid = 20, nvgrid = 20)
#' \donttest{
#' # complete-case estimation discards rows with a missing mark; also, no
#' auxiliary covariate is needed
#' fit <- kernel_sievePH(eventTime, eventInd, mark, tx, nboot = 50,
#'                           missmethod = "CC", tau = 3, tband = 0.5, hband = 0.3,
#'                           a = 0.1, b = 1, ntgrid = 20, nvgrid = 20)
#' }
#'
#' @importFrom plyr laply
#'
#' @export
kernel_sievePH <- function(eventTime, eventInd,mark, tx, aux = NULL , strata = NULL, nboot = 500, missmethod,
                               tau, tband, hband, a, b, ntgrid = 100, nvgrid = 100, confLevel = 0.95, seed = 62648) {



  # nauxiliary a numeric indicator for whether any auxiliary variable is used in augmenting the IPW estimation equations
  if(missmethod == "CC"){
    nauxiliary = 0
    missmark = 0
  }else if (missmethod == "IPW"){
    nauxiliary = 0
    missmark = 1
  }else if (missmethod == "AIPW"){
    nauxiliary = 1
    missmark = 1
  }else{
    print("Missing-Mark Method Missing")
  }

  a0 <- a
  b1 <- b
  nsamp <- NULL
  # kk is the number of strata
  if(is.null(strata)){
    nsamp[1] = length(eventTime)
    kk <- 1
  }else{
    stratum <- unique(strata)
    kk <- length(stratum)
    for(i in 1:kk){
      nsamp[i] <- sum(strata == stratum[i])
    }
  }

  ncov = 1 #number of covariates
  ekconst <- 3/5

  set.seed(seed)

  #*****************************************************************************
  #     Estimation and Hypothesis tests for Mark PH model with missing marks
  #
  #*****************************************************************************



 # IMPORTANT NOTE

 # This version of the code inputs mark on its own scale, and re-scales it to 0-1 before doing the analysis.
  vV <- mark
  mn <- min(mark, na.rm = TRUE) #minimum observed mark value
  mx <- max(mark, na.rm = TRUE) #maximum observed mark value
  # Put the marks on the scale 0 to 1:
  vV <- (mark - mn)/(mx-mn)
  # Set mark with missing values to be a large number
  vV[is.na(vV)] <- 8888

  rinfect1 <- sum(eventInd>=0.5 & tx >=0.5)
  rinfect0 <- sum(eventInd>=0.5 & tx < 0.5)

  #total sample
  nsampa <- sum(nsamp)
  rnsamp = nsampa
  #*****************************************************************************
  #                      Set up the smoothing parameters
  #
  #*****************************************************************************

  # Set up the smoothing parameters
  tstep <- tau / as.numeric(ntgrid)
  vstep <- 1 / as.numeric(nvgrid)
  iskip <- a0 / vstep

  # Epanechnikov kernel has squared integral = 3/5 !
  ekconst <- 3/5

  # calculate the test statistics for H_{10}:VE(v)=0 for v over [a1,b1], 0<a1<b1=<1
  # calculate the test statistics for H_{20}:VE(v)=VE for v over [a2,b1], a2>a1 >=a0

  a1 <- a0 + vstep
  # a1 = a' in Gilbert and Sun (2014, figure 3 caption)
  a2 <- a1 + 2 * vstep


  #*****************************************************************************
  #                     Set up the control parameters
  #
  #     For the procedures without missing marks, set MissMark=0;
  #     For the procedures with missing marks, set MissMark=1.
  #      (MissMark=0) corresponds to setting R(ks,i)=1 and wght(ks,i)=1.
  #
  #     For the procedure without using the auxiliary marks, set Nauxiliary=0;
  #        to use the auxiliary marks, set Nauxiliary=1.
  #*****************************************************************************


  # Read data into the designed variables

  sumdelta <- 0
  Rsum <- 0

  #covart does not allow time-varying covariates in this version
  covart <- array(0, dim = c(kk, ncov, max(nsamp)))
  time <- matrix(0, nrow = kk, ncol = max(nsamp))
  censor <- matrix(0, nrow = kk, ncol = max(nsamp))
  markm <- matrix(0, nrow = kk, ncol = max(nsamp))
  Ax <- matrix(0, nrow = kk, ncol = max(nsamp))
  R <- matrix(0, nrow = kk, ncol = max(nsamp))

 if(kk == 1){
    n_ks <- nsamp[1]
    ks = 1
    covart[ks, 1, 1:n_ks] <- tx[1:n_ks]
    time[ks, 1:n_ks] <- eventTime[1:n_ks]
    censor[ks, 1:n_ks] <- eventInd[1:n_ks]
    Ax[ks, 1:n_ks] <- aux[1:n_ks]
    R[ks, 1:n_ks] <- 1-is.na(mark[1:n_ks])
    # Important For all eventInd ==0, R is set as 1
    R[ks,eventInd==0] <- 1
    markm[ks, 1:n_ks] <- vV[1:n_ks]


  }else{
    for (ks in 1:kk) {
      n_ks <- nsamp[ks]
      indice_ks <- (strata == stratum[ks])
      covart[ks, 1, 1:n_ks] <- tx[indice_ks]
      time[ks, 1:n_ks] <- eventTime[indice_ks]
      censor[ks, 1:n_ks] <- eventInd[indice_ks]
      markm[ks, 1:n_ks] <- vV[indice_ks]
      Ax[ks, 1:n_ks] <- aux[indice_ks]
      R[ks, 1:n_ks] <- 1-is.na(mark[indice_ks])
      # Important For all eventInd ==0, R is set as 1
      R[ks,censor[ks, 1:n_ks]==0] <- 1

    }
  }




  # Calculate sumdelta and Rsum
  sumdelta <- sum(censor)
  Rsum <- sum(R * censor)

  # Calculate the probability of nonmissing marks.

  # wght(ks,i) can be calculated using an external program

  # setup the covariates for logit regression for R:

  # Assuming kk, nsamp, time, censor, and missmark are already defined

  # Define ncovR
  ncovR <- 2

  # Create covartR array
  n_max <- max(nsamp)
  covartR <- array(0, dim = c(kk, ncovR, n_max))
  for(ks in 1:kk){
    covartR[ks,1,1:nsamp[ks]] <- 1
    covartR[ks,2,1:nsamp[ks]] <- time[ks,1:nsamp[ks]]

  }

  psiRse <- matrix(0, nrow = kk, ncol = ncovR)
  # Handle missmark
  if (missmark == 0) {
    R <- matrix(1, nrow = kk, ncol = n_max)
    wght <- R
  } else {

    estplogitans <- estplogit(kk, nsamp, ncovR, covartR, R, censor)
    psiR <- estplogitans$PSI
    psiRvar <- estplogitans$FVAR

    for(ks in 1:kk){
      psiRse[ks,] <- sqrt(diag(psiRvar[ks,,]))
    }


    #logistic regression
    # xx = R[1,censor[1,]==1]
    # fit = summary(glm(summary(glm(xx~covartR[1,2,censor[1,]==1]),family = binomial(link = "logit"))))
    wght <- matrix(0, nrow = kk, ncol = n_max)
    for(ks in 1:kk){
     ei <- covartR[ks, 1, ] * psiR[ks, 1] + covartR[ks, 2, ] * psiR[ks, 2]
     phi <- exp(ei) / (1 + exp(ei))
     wght[ks,] <- R[ks,] / (censor[ks,] * phi + 1 - censor[ks,])
   }


  }


  # the end of calculating wght(ks,i)

  # Calculate the pdf G(ks,i,ispot) for the second stage estimator, aipw

  # IF use (Nauxiliary.eq.1),
  # the pdf G(ks,i,ispot) can also be calculated using an external program

  # This version of the code only handles a special case where the auxiliary for
  # association with the mark variable V is binary; and the failure time is the variable for
  # predicting whether a case has the mark V available.  These parts of the program would
  # need to be changed manually for other types of auxiliary covariates.
  # The relevant code is surrounded by the comments AUXILIARY CODE

  # AUXILIARY CODE
  # gender is the only significant variable for the mark.
  # We set gender as the only auxliary variable.
  # Since gender is a binary variable, do the logistic regression
  # of A(gender) ~ V,T,Z .

  # In addition, the time is the only significant variable for predicting
  # the probability of missingness.
  # We use the logistic regression of R ~ T.

  # setup the covariates for logistic regression for gender to calculate g(a|w):


  # Define ncovG
  ncovG <- 4

  # Get maximum sample size across all strata
  n_max <- max(nsamp)

  # Create covartG array
  covartG <- array(0, dim = c(kk, ncovG, n_max))
  for(ks in 1:kk){
    covartG[ks,1,1:nsamp[ks]] <- 1
    covartG[ks,2,1:nsamp[ks]] <- time[ks,1:nsamp[ks]]
    covartG[ks,3,1:nsamp[ks]] <- covart[ks,1,1:nsamp[ks]]#treatment
    covartG[ks,4,1:nsamp[ks]] <- markm[ks,1:nsamp[ks]]
  }

  # Compute cenG
  cenG <- censor * R

  if (nauxiliary == 1) {
    G <- array(1, dim = c(kk, n_max, nvgrid))

    #estplogit(kk, nsamp, ncovG, covartG, Ax, cenG, psiG, psiGvar)
    estplogitGans <- estplogit(kk, nsamp, ncovG, covartG, Ax, cenG)
    psiG <- estplogitGans$PSI
    psiGvar <- estplogitGans$FVAR

    psiGse <- matrix(0, nrow = kk, ncol = ncovG)

    for(ks in 1:kk){
      psiGse[ks,] <- sqrt(diag(psiGvar[ks,,]))
    }

    for(ks in 1:kk){
      for(i in 1:nsamp[ks]){
        ei = 0
        for(j in 1:(ncovG-1)){
          ei <- ei + covartG[ks,j,i]*psiG[ks,j]
        }

        vspot <- (1:nvgrid)*vstep
        phi <- exp(ei + psiG[ks, ncovG] * vspot) / (1+exp(ei + psiG[ks, ncovG] * vspot))
        G [ks,i,1:nvgrid]<- phi * Ax[ks,i] + (1 - phi) * (1 - Ax[ks,i])

      }
    }
  } else {
    G <- array(1, dim = c(kk, n_max, nvgrid))
  }

  #browser()
  # the end of calculating G(ks,i,ispot).
  # AUXILIARY CODE

  # Estimation of the beta(v) at each mark v

  # Comparison with standard PH model analysis
  # Also serve as an initial value for the mark-specific PH model

  # CALL COX function


  # library(survival)
  # phfit <- summary(coxph(Surv(time[1,], censor[1,])~trt))
  coxfit <- estpvry (tau, kk, nsamp, ncov, time, covart, censor)
  betacox <- coxfit$BETA
  varcox <- coxfit$var


  # Estimate of s.e. of beta:
  seph <- sqrt(varcox[1, 1])

  # Wald test for beta=0:
  Tw <- betacox[1] / seph

  # uniker is (1/h)K((t_k-t)/h) delta where K is uniform kernel
  # normker is (1/h)K((t_k-t)/h) delta where K is normal kernel

  betacom <- matrix(0,ncol = nvgrid,nrow = ncov)
  secom <- matrix(0,ncol = nvgrid,nrow = ncov)

  betaipw <- matrix(0,ncol = nvgrid,nrow = ncov)
  seipw <- matrix(0,ncol = nvgrid,nrow = ncov)
  LAMBDA0 <- array(0, dim = c(kk, max(nsamp), nvgrid))

  betaaug <- matrix(0,ncol = nvgrid,nrow = ncov)
  seaug <- matrix(0,ncol = nvgrid,nrow = ncov)

  betaofv <- array(0, dim = c(ncov, nvgrid))
  SigmaInv <- array(0, dim = c(ncov, ncov,nvgrid))

  # baseline hazard functions:
  Lamktv0 <- array(0, dim = c(kk, ntgrid, nvgrid))

  #VE
  veest <- NULL

  vese <- NULL


  ###################################################################
  for (ispot in 1:nvgrid) {
    vspot <- ispot * vstep
    deltak <- matrix(0,nrow = kk, ncol = max(nsamp))
    #function(tk, tvalue, hband, delt)
    for(ks in 1:kk){
      deltak[ks,] <- Epanker(markm[ks,], vspot, hband, censor[ks,])
    }


    # Complete data likelihood estimator: complete
    # Set the weight wght(ks, i) = R(ks, i) to get the complete data likelihood estimator
    #browser()
    #comfit <- estpcom(tau, kk, nsamp, ncov, time, covart, deltak, R,betacox)
    comfit <- estpcomcplusplus(tau, kk, nsamp, ncov, time, covart, deltak, R,betacox)

    betacom[,ispot] <- comfit$BETA
    secom[,ispot] <- sqrt(diag(comfit$var))


    # The inverse probability weighted estimator: ipw

    #ipwfit <- estpipw(tau, kk, nsamp, tband, ncov, time, covart, censor, deltak, wght, betacom[,ispot])
    ipwfit <- estpipwcplusplus(tau, kk, nsamp, tband, ncov, time, covart, censor, deltak, wght, betacom[,ispot])

    betaipw[,ispot] <- ipwfit$BETA
    seipw[,ispot] <- sqrt(diag(ipwfit$var))
    LAMBDA0[,,ispot] <- ipwfit$LAMBDA0
  }



  #   ********************************************************************
  #                         second stage estimator:   aipw
  #   ********************************************************************

    # Second stage estimator: aipw
    CLAMBDAIPW <- array(0, dim = c(kk, max(nsamp)))
    LAMBDAIPW <- array(0, dim = c(kk, max(nsamp), nvgrid))

    LAMBDAUG<- array(0, dim = c(kk, max(nsamp), nvgrid))
    CLAMBDAUG<- array(0, dim = c(kk, max(nsamp)))

    for (ispot in 1:nvgrid) {
      vspot <- ispot * vstep
      # CALCULATE CLAMBDA(KS,I)= $\int_0^1 \lambda_k(t,u|z) h_k(a|t,u,z) du$ AT $W_{ki}=(T_{ki},Z_{ki})$:
      for (ks in 1:kk) {
        valid_indices <- which((time[ks, 1:nsamp[ks]] <= tau) & (censor[ks,1:nsamp[ks] ] > 0.5))

        #vectorize over loop I = 1, NSAMP(ks)
        if (any(valid_indices)) {

          if(ncov==1){
            BZIPW <- NULL
            BZIPW[valid_indices] <- covart[ks,1 , valid_indices] * betaipw[1, ispot]
            BZIPW <- as.matrix(BZIPW, ncol = 1)
          }else{
            BZIPW <- array(0, dim = c(1, max(kk)))
            BZIPW[valid_indices,1] <- t(covart[ks, , valid_indices]) %*% betaipw[, ispot] #summing over J

          }

          #the original fortran code had LAMBDA0 for nvgrid only
          LAMBDAIPW[ks, valid_indices, ispot] <- LAMBDA0[ks, valid_indices,ispot] * exp(BZIPW[valid_indices,1])
          CLAMBDAIPW[ks, valid_indices] <- CLAMBDAIPW[ks, valid_indices] + LAMBDA0[ks, valid_indices,ispot] * exp(BZIPW[valid_indices,1]) * G[ks,valid_indices,ispot]*vstep
        }
      }
    }
#Note: LAMBDAIPW is lambda_k(t,u|z) in equation 11, Sun and Gilbert 2012
#      CLAMBDAIPW is the int_0^1 lambda_k(t,u|z)g_k(a|t,u,z)du in equation 11
########################################################
# CALCULATE CKLAMBDA(KS,I, u)= $ K_h(u-v) \lambda_k(t,u|z) h_k(a|t,u,z) du$ AT $W_{ki}=(T_{ki},Z_{ki})$:
# K_h(u-v)du is
#DRHOipw is the derivative of rho_k^{ipw}(v,w_i) w.r.t v in equation 12. I don't understand why kernel-matrix needs to be in the calculation for CKLAMBDAIPW
    CKLAMBDAIPW <- array(0, dim = c(kk, max(nsamp),nvgrid))
    for (ks in 1:kk) {
      for (i in 1:nsamp[ks]) {
        vspot <- (1:nvgrid) * vstep

        kernel_matrix <- laply(vspot, function(x){
          z <- EpankerV(vspot, x, hband, censor[ks, i])
        })

        CKLAMBDAIPW[ks, i, ] <- kernel_matrix %*% (G[ks,i,]*LAMBDAIPW[ks, i, ]) * vstep
      }
    }

    CLAMBDAUG <- matrix(0, nrow = kk, ncol = max(nsamp))

    CUMB1 <- rep(0,nvgrid)

    vspot <- (1:nvgrid) * vstep

    for (ispot in 1:nvgrid) {

      vspot <- ispot * vstep
      DRHOipw <- matrix(0,ncol = max(nsamp),nrow = kk)
      deltak <- matrix(0,ncol = max(nsamp),nrow = kk)
      for (ks in 1:kk) {
        valid_indices <- which((time[ks, 1:nsamp[ks]] <= tau) & (censor[ks, 1:nsamp[ks]] > 0.5) & (CLAMBDAIPW[ks, 1:nsamp[ks]] > 0.00000001))
        deltak[ks, ] <- Epanker(markm[ks, ], vspot, hband, censor[ks, ])
        DRHOipw[ks, valid_indices] <- CKLAMBDAIPW[ks, valid_indices, ispot] / CLAMBDAIPW[ks, valid_indices]
        DRHOipw[ks, !valid_indices] <- 0.0
      }

      # initial value
      betaaug0 <- betaipw[, ispot]


      #estpaug_result <- estpaug(tau, tstep, ntgrid, tband, kk, nsamp, ncov, time, covart, censor, deltak, wght, DRHOipw, betaaug0)
      estpaug_result = estpaugcplusplus(tau, tstep, ntgrid, tband, kk, nsamp, ncov, time, covart, censor, deltak, wght, DRHOipw, betaaug0)
      betaaug[, ispot] <- estpaug_result$BETA
      seaug[, ispot] <- sqrt(diag(estpaug_result$var))



      if (ispot >= iskip & ispot > 1) {
        # calculate cumulative beta_1(v):
        CUMB1[ispot] <- CUMB1[ispot - 1] + betaaug[1, ispot] * vstep
      }else if(ispot >= iskip & ispot == 1){
        CUMB1[ispot] <- betaaug[1, ispot] * vstep
      }

      # Inverse matrix of information matrix:
      betaofv[, ispot] <- betaaug[, ispot]
      SigmaInv[, , ispot] <- estpaug_result$F

      # calculate the baseline hazard functions:
      Lamktv0[, , ispot] <- estpaug_result$LAMBDAk0
      LAMBDA0 <- estpaug_result$LAMBDA0
      # CALCULATE CLAMBDA(KS,I)= $\int_0^1 \lambda_k(t,u|z) h_k(a|t,u,z) du$ AT $W_{ki}=(T_{ki},Z_{ki})$:
      for (ks in 1:kk) {
        valid_indices <- which((time[ks, 1:nsamp[ks]] <= tau) & (censor[ks, 1:nsamp[ks]] > 0.5))

        if(ncov==1){
          BZaug <- as.matrix(covart[ks, , valid_indices] * betaaug[1, ispot],ncol=1)
        }else{
          BZaug <- covart[ks, , valid_indices] %*% betaaug[, ispot]
        }


        LAMBDAUG[ks, valid_indices, ispot] <- LAMBDA0[ks, valid_indices] * exp(BZaug)
        CLAMBDAUG[ks, valid_indices] <- CLAMBDAUG[ks, valid_indices] + G[ks, valid_indices, ispot] * LAMBDA0[ks, valid_indices] * exp(BZaug) * vstep
      }

      # estimate of VE(v):
      veest[ispot] <- 1 - exp(betaaug[1, ispot])

      # One estimate of s.e. of \hat VE(v):
      vese[ispot] <- seaug[1, ispot] * exp(betaaug[1, ispot])


    }
    # 95% confidence interval for beta1(v):
    cv_confLevel <- qnorm(1-(1-confLevel)/2)
    beta1lw=betaaug[1,]-cv_confLevel*seaug[1,]
    beta1up=betaaug[1,]+cv_confLevel*seaug[1,]
    markgrid <- (1:nvgrid) * vstep
    velw <- 1-exp(beta1up)
    veup <- 1-exp(beta1lw)
    te <- data.frame("mark" = markgrid*(mx-mn)+mn, "TE" = veest,"SE" = vese, "LB" = velw, "UB" = veup)




    # **************************************************************************
    # calculate the quantities needed for GDIST2N
    # **************************************************************************
    AsigInv <- array(0, dim = c(ncov*ncov, nvgrid, nvgrid))
    tempaug <- array(0, dim = c(kk, max(nsamp), nvgrid))
    for (ispot in iskip:nvgrid) {
      #*** re: delete a+
      # va=a0+vstep*(float(ispot)-1.0)
      # vb=a0+vstep*float(ispot)
      # vspot=a0+float(ispot)*vstep

      va <- vstep * (ispot - 1.0)
      vb <- vstep * ispot
      vspot <- ispot * vstep

      sigtemp <- matrix(0.0, nrow = ncov, ncol = ncov)
      for(mspot in iskip:nvgrid){
        vmspot <- mspot * vstep
        sigtemp <- sigtemp + SigmaInv[, , mspot] * Epanker(vspot, vmspot, hband, 1.0) * vstep
        AsigInv[, ispot, mspot] <- as.vector(sigtemp[1:ncov,1:ncov])
      }



      DRHOaug <- matrix(0,nrow = kk, ncol = max(nsamp))
      for(ks in 1:kk){
        valid_indices <- which((time[ks,] <= tau) & (censor[ks,] > 0.5) & (CLAMBDAUG[ks,] > 0.00000001))
        DRHOaug[ks,valid_indices] <- G[ks,valid_indices, ispot] * LAMBDAUG[ks,valid_indices, ispot] / CLAMBDAUG[ks,valid_indices]
        DRHOaug[ks,!valid_indices] <- 0.0

        temp1 <- (1.0 - wght[ks,]) * DRHOaug[ks,] * vstep
        temp1[!(time[ks,] <= tau) | !(censor[ks,] > 0.5)] <- 0.0

        temp2 <- wght[ks,]
        temp2[!(time[ks,] <= tau) | !(markm[ks,] > va) | !(markm[ks,] <= vb) | !(censor[ks,] > 0.5)] <- 0.0

        tempaug[ks,1:nsamp[kk], ispot] <- temp1 + temp2
      }



    }

#############################################
    S0N <- array(0,dim = c(kk, max(nsamp),nvgrid))
    S1N <- array(0,dim = c(kk*ncov, max(nsamp), nvgrid))

    for (ispot in iskip:nvgrid) {
      for (ks in 1:kk) {
        # for time independent covariates
        BZt <- matrix(0, nrow = nsamp[ks], ncol = 1)

        if(ncov==1){
          BZt[1:nsamp[ks],1] <- as.matrix(covart[ks, , 1:nsamp[ks]] * betaofv[1, ispot], ncol = 1)

        }else{
          BZt <- t(covart[ks, , 1:nsamp[ks]]) %*% betaofv[, ispot]
        }


        S0 <- rep(0, nsamp[ks])
        S1 <- matrix(0, nrow = ncov, ncol = nsamp[ks])

        for (i in 1:nsamp[ks]) {
          valid_indices <- which(time[ks, 1:nsamp[ks]] >= time[ks, i])
          S0[i] <- sum(exp(BZt[valid_indices,1]))
          for(J in 1:ncov){
            S1[J,i] <- sum(covart[ks,J , valid_indices]*exp(BZt[valid_indices,1]))
          }
        }


        S0N[ks, 1:nsamp[ks], ispot] <- S0
        arrayseq <- ((ks-1)*ncov +1): ks*ncov
        S1N[ arrayseq,1:nsamp[ks], ispot] <- S1
      }
    }


    # **************************************************************************
    # finish calculating the quantities needed for GDIST2N
    # **************************************************************************

    # **************************************************************************
    # Simulate distribution of \sqrt n(\hat B(v)-B(v)) based on aipw
    # **************************************************************************



    CUMB1x <- rep(0, nvgrid)
    CUMB1se <- rep(0, nvgrid)
    BootDist <- matrix(0,nrow = nboot, ncol = nvgrid)
    for (iboot in 1:nboot) {
      zdev <- matrix(rnorm(kk*max(nsamp)), nrow = kk, ncol = max(nsamp))

      #browser()
      #CUMBDIST is simulated B(v)
      CUMBDIST <- GDIST2Ncplusplus(nvgrid, iskip, zdev, kk, nsamp, ncov, time, covart,
                                    betaofv, SigmaInv, S0N, S1N, tempaug, AsigInv)


      # calculate the Gaussian multiplier process CUMBDIST for each iboot
      #CUMBDIST <- GDIST2N(nvgrid, iskip, zdev, kk, nsamp, ncov, time, covart,
      #                    betaofv, SigmaInv, S0N, S1N, tempaug, AsigInv)

      CUMB1x[iskip:nvgrid] <- CUMB1x[iskip:nvgrid] + CUMBDIST[1, iskip:nvgrid]
      CUMB1se[iskip:nvgrid] <- CUMB1se[iskip:nvgrid] + CUMBDIST[1, iskip:nvgrid]^2.0

      BootDist[iboot, ] <- CUMBDIST[1, ]
    }



    CUMB1x[iskip:nvgrid] <- CUMB1x[iskip:nvgrid] / nboot
    CUMB1se[iskip:nvgrid] <- sqrt((CUMB1se[iskip:nvgrid] - (CUMB1x[iskip:nvgrid])^2.0 * nboot) / nboot)
    # do ispot=iskip,nvgrid
    # for (ispot in iskip:nvgrid) {
    #   # bootstrap standard error of the IPWA estimator $\hat B_1(v)-B_1(v)$:
    #   CUMB1x[ispot] <- CUMB1x[ispot] / nboot
    #   CUMB1se[ispot] <- sqrt((CUMB1se[ispot] - (CUMB1x[ispot])^2.0 * nboot) / nboot)
    # }


    # **************************************************************
    # This data application version is different from the simulation
    # version; it is more flexible; it can construct test statistics
    # on [a1,b1] or [a2,b1].
    # **************************************************************

    # calculate the test statistics for H_{10}:VE(v)=0 for v over [a1,b1], b1=1
    iskipa1 <- a1 / vstep
    nvgrid1 <- b1 / vstep

    # calculate the test statistics for H_{20}:VE(v)=VE for v over [a2,b1], a2>a1
    iskipa2 <- a2 / vstep

    cBproc1 <- rep(0, nvgrid)
    cBproc2 <- rep(0, nvgrid)

    cBproc1[iskipa1] <- CUMB1[iskipa1] - CUMB1[iskipa1]
    cBproc2[iskipa2] <- (CUMB1[iskipa2] - CUMB1[iskipa1]) / (a2 - a1) - (CUMB1[nvgrid1] - CUMB1[iskipa1]) / (b1 - a1)

    TSUP1 <- abs(cBproc1[iskipa1])
    TSUP1m <- cBproc1[iskipa1]
    TSUP2 <- abs(cBproc2[iskipa2])
    TSUP2m <- cBproc2[iskipa2]
    Tint1 <- 0
    Tint1m <- 0
    Tint2 <- 0
    Tint2m <- 0
    #browser()
    for (ispot in (iskipa1 + 1):nvgrid1) {
      vspot <- ispot * vstep

      cBproc1[ispot] <- CUMB1[ispot] - CUMB1[iskipa1]
      cBproc2[ispot] <- (CUMB1[ispot] - CUMB1[iskipa1]) / (vspot - a1) - (CUMB1[nvgrid1] - CUMB1[iskipa1]) / (b1 - a1)

      # two-sided tests for testing $H_1$:
      TSUP1 <- max(TSUP1, abs(cBproc1[ispot]))
      Tint1 <- Tint1 + cBproc1[ispot]^2 * (CUMB1se[ispot]^2 - CUMB1se[ispot - 1]^2)

      # one-sided tests for testing $H_1$:
      TSUP1m <- min(TSUP1m, cBproc1[ispot])
      Tint1m <- Tint1m + cBproc1[ispot] * (CUMB1se[ispot]^2 - CUMB1se[ispot - 1]^2)

      if (ispot >= iskipa2) {
        # two-sided tests for testing $H_2$:
        TSUP2 <- max(TSUP2, abs(cBproc2[ispot]))
        Tint2 <- Tint2 + cBproc2[ispot]^2 * (CUMB1se[ispot]^2 - CUMB1se[ispot - 1]^2)

        # one-sided tests for testing $H_2$:
        TSUP2m <- min(TSUP2m, cBproc2[ispot])
        Tint2m <- Tint2m + cBproc2[ispot] * (CUMB1se[ispot]^2 - CUMB1se[ispot - 1]^2)
      }
    }


    # *****************************************************************
    # Find the p-values of the test statistics.
    # *****************************************************************

    # the initial values for p-values
    TSUP1pv <- 0.0
    TSUP1mpv <- 0.0
    TSUP2pv <- 0.0
    TSUP2mpv <- 0.0
    Tint1pv <- 0.0
    Tint1mpv <- 0.0
    Tint2pv <- 0.0
    Tint2mpv <- 0.0
    cBproc1s <- matrix(0,nrow = nvgrid, ncol = nboot)
    cBproc2s <- matrix(0,nrow = nvgrid, ncol = nboot)

    for (iboot in 1:nboot) {
      cBproc1x <- rep(0, nvgrid)
      cBproc2x <- rep(0, nvgrid)

      cBproc1x[iskipa1] <- BootDist[iboot, iskipa1] - BootDist[iboot, iskipa1]
      temp2nd <- (BootDist[iboot, nvgrid1] - BootDist[iboot, iskipa1]) / (b1 - a1)
      cBproc2x[iskipa2] <- (BootDist[iboot, iskipa2] - BootDist[iboot, iskipa1]) / (a2 - a1) - temp2nd

      TSUP1x <- abs(cBproc1x[iskipa1])
      TSUP1mx <- cBproc1x[iskipa1]
      TSUP2x <- abs(cBproc2x[iskipa2])
      TSUP2mx <- cBproc2x[iskipa2]
      Tint1x <- 0
      Tint1mx <- 0
      Tint2x <- 0
      Tint2mx <- 0

      for (ispot in (iskipa1 + 1):nvgrid1) {
        vspot <- ispot * vstep

        cBproc1x[ispot] <- BootDist[iboot, ispot] - BootDist[iboot, iskipa1]
        cBproc2x[ispot] <- (BootDist[iboot, ispot] - BootDist[iboot, iskipa1]) / (vspot - a1) - temp2nd

        # for plotting the Gaussian multiplier realizations; the actual processes need to multiply by sqrt(n):
        cBproc1s[ispot, iboot] <- sqrt(rnsamp) * cBproc1x[ispot]
        cBproc2s[ispot, iboot] <- sqrt(rnsamp) * cBproc2x[ispot]

       # two-sided tests for testing $H_1$:
        TSUP1x <- max(TSUP1x, abs(cBproc1x[ispot]))
        Tint1x <- Tint1x + cBproc1x[ispot]^2 * (CUMB1se[ispot]^2 - CUMB1se[ispot - 1]^2)

        # one-sided tests for testing $H_1$:
        TSUP1mx <- min(TSUP1mx, cBproc1x[ispot])
        Tint1mx <- Tint1mx + cBproc1x[ispot] * (CUMB1se[ispot]^2 - CUMB1se[ispot - 1]^2)

        if (ispot >= iskipa2) {
          # two-sided tests for testing $H_2$:
          TSUP2x <- max(TSUP2x, abs(cBproc2x[ispot]))
          Tint2x <- Tint2x + cBproc2x[ispot]^2 * (CUMB1se[ispot]^2 - CUMB1se[ispot - 1]^2)

          # one-sided tests for testing $H_2$:
          TSUP2mx <- min(TSUP2mx, cBproc2x[ispot])
          Tint2mx <- Tint2mx + cBproc2x[ispot] * (CUMB1se[ispot]^2 - CUMB1se[ispot - 1]^2)
        }
      }

      # calculate bootstrap p-values for H10:
      TSUP1pv <- TSUP1pv + as.numeric(TSUP1x > TSUP1) / nboot
      TSUP1mpv <- TSUP1mpv + as.numeric(TSUP1mx < TSUP1m) / nboot
      Tint1pv <- Tint1pv + as.numeric(Tint1x > Tint1) / nboot
      Tint1mpv <- Tint1mpv + as.numeric(Tint1mx < Tint1m) / nboot

      # calculate bootstrap p-values for H20:
      TSUP2pv <- TSUP2pv + as.numeric(TSUP2x > TSUP2) / nboot
      TSUP2mpv <- TSUP2mpv + as.numeric(TSUP2mx < TSUP2m) / nboot
      Tint2pv <- Tint2pv + as.numeric(Tint2x > Tint2) / nboot
      Tint2mpv <- Tint2mpv + as.numeric(Tint2mx < Tint2m) / nboot


    }


    # the actual test statistics need to multiply by sqrt(n):
    TSUP1 <- sqrt(rnsamp) * TSUP1
    TSUP1m <- sqrt(rnsamp) * TSUP1m
    Tint1 <- sqrt(rnsamp) * Tint1
    Tint1m <- sqrt(rnsamp) * Tint1m
    TSUP2 <- sqrt(rnsamp) * TSUP2
    TSUP2m <- sqrt(rnsamp) * TSUP2m
    Tint2 <- sqrt(rnsamp) * Tint2
    Tint2m <- sqrt(rnsamp) * Tint2m


    # for plotting the test processes and the Gaussian multiplier processes;
    # the actual test processes need to multiply by sqrt(n):
    for (ispot in (iskipa1 + 1):nvgrid) {
      vspot <- ispot * vstep
      cBproc1[ispot] <- sqrt(rnsamp) * cBproc1[ispot]

      if (ispot >= iskipa2) {
        cBproc2[ispot] <- sqrt(rnsamp) * cBproc2[ispot]

      }
    }

    Rmiss <- Rsum / sumdelta
    #browser()
    test10 <- c(TSUP1,TSUP1m, Tint1,Tint1m)
    pvalue10 <- c(TSUP1pv,TSUP1mpv, Tint1pv,Tint1mpv)
    ans10 <- data.frame(rbind(test10,pvalue10))
    colnames(ans10) <- c("TSUP1","TSUP1m", "Tint1","Tint1m")
    rownames(ans10) <- c("Test Statistic", "P-value")

    test20 <- c(TSUP2,TSUP2m, Tint2,Tint2m)
    pvalue20 <- c(TSUP2pv,TSUP2mpv, Tint2pv,Tint2mpv)
    ans20 <- data.frame(rbind(test20,pvalue20))
    colnames(ans20) <- c("TSUP2","TSUP2m", "Tint2","Tint2m")
    rownames(ans20) <- c("Test Statistic", "P-value")


    cBproc1.df <- data.frame(cbind(te$mark[(iskipa1 + 1):nvgrid],vstep*((iskipa1 + 1):nvgrid),cBproc1[(iskipa1 + 1):nvgrid],cBproc1s[(iskipa1 + 1):nvgrid,]))
    colnames(cBproc1.df) <- c("Mark", "Standardized mark","Observed", paste0("S", 1:nboot))
    cBproc2.df <- data.frame(cbind(te$mark[(iskipa2):nvgrid],vstep*((iskipa2):nvgrid),cBproc2[iskipa2:nvgrid], cBproc2s[iskipa2:nvgrid, ]))
    colnames(cBproc2.df) <- c("Mark", "Standardized mark", "Observed", paste0("S", 1:nboot))
    return(list("cBproc1" = cBproc1.df,
                "cBproc2" = cBproc2.df,
                "H10" = ans10,
                "H20" = ans20,
                "te" = te,
                "estBeta" = data.frame("mark" = te$mark,"betacom" = t(betacom),
                                   "secom" = t(secom),
                                   "betaipw" = t(betaipw),
                                   "seipw" = t(seipw),
                                   "betaaug" = t(betaaug),
                                   "seaug" = t(seaug))
                ))

}


# **************************************************************************
#
# function matrix A times B gives C
# NP, MP, LP are the declared dimensions and N, M, L are
# the actual dimensions
# **************************************************************************

ATIMEB <- function(A, B, N, M, L) {
  C <- A[1:N, 1:M] %*% B[1:M, 1:L]
  return(C)
}

# C     *****************************************************************
# C
# C                         _ESTPVRY       function_
# C
# C     _estimate the time-varying coefficient in the cox model_
# C     *****************************************************************

estpvry <- function(tau, KK, N, NP, X, ZT, DELTA) {

  BETA = rep(0, NP)
  BETA[1] = -0.7
  var <- matrix(0, nrow = NP, ncol = NP)

  KL <- 6

  for (KITER in 1:KL) {
    U <- numeric(NP)
    F <- matrix(0, nrow = NP, ncol = NP)
    F2 <- matrix(0, nrow = NP, ncol = NP)

    for (ks in 1:KK) {
      if(NP==1){
        BZt <- as.matrix(ZT[ks, , 1:N[ks]] * BETA, ncol = 1)

      }else{
        BZt <- ZT[ks, , 1:N[ks]] %*% BETA

      }

      # c*******************************************************************

      mask <- X[ks, 1:N[ks]] <= tau
      indices <- which(mask)

      for (i in indices) {
        S0 <- numeric(N[ks])
        S1 <- matrix(0, nrow = NP, ncol = N[ks])
        S2 <- matrix(0, nrow = NP, ncol = NP)

        L <- which(X[ks, 1:N[ks]] >= X[ks, i])


        if(NP == 1){
          ZTL <- ZT[ks,1:NP , L]
          S0[i] <- sum(exp(BZt[L,1]))
          S1[1, i] <- sum(ZTL*exp(BZt[L,1]))
          S2 <- S2 + sum(ZTL*ZTL*exp(BZt[L,1]))

        }else{
          S0[i] <- sum(exp(BZt[L, 1]))
          S1[1:NP, i] <-  Zt[ks, 1:NP, L] %*% exp(BZt[L, 1])
          for(J in 1:NP){
            for(K in 1:NP){
              S2[J, K] <- t(exp(BZt[L, 1]))%*%(Zt[ks,J , L]*Zt[ks,K , L])
            }
          }
        }

        # c???
        if (S0[i] > 0.00000001) {
          U <- U + DELTA[ks, i] * (ZT[ks, , i] - S1[, i] / S0[i])
          F <- F + DELTA[ks, i] * (S2 / S0[i] - tcrossprod(S1[, i]) / S0[i]^2)
          F2 <- F2 + DELTA[ks, i]^2 * (S2 / S0[i] - tcrossprod(S1[, i]) / S0[i]^2)
        }
      }
    }

    BETA <- BETA + solve(F, U)
  }

  # C     the variance of \hat\beta(v) is var in the following
  var <- solve(F) %*% F2 %*% solve(F)

  return(list(BETA = BETA, var = var))
}

# C     *****************************************************************
# C
# C                         _ESTPCOM       function_
# C
# C     _Complete data estimate the varying coefficient in the mark PH model_
# C     *****************************************************************
estpcom <- function(tau, KK, N, NP, X, ZT, DELTA, WGHT,beta0) {
  #beta0 is initial value
  BETA = beta0

  var <- matrix(0, nrow = NP, ncol = NP)

  KL <- 6

  for (KITER in 1:KL) {
    U <- numeric(NP)
    F <- matrix(0, nrow = NP, ncol = NP)
    F2 <- matrix(0, nrow = NP, ncol = NP)

    for (ks in 1:KK) {

      if(NP==1){
        BZt <- as.matrix(ZT[ks, , 1:N[ks]] * BETA, ncol = 1)

      }else{
        BZt <- ZT[ks, , 1:N[ks]] %*% BETA}


      mask <- X[ks, 1:N[ks]] <= tau
      indices <- which(mask)
      S0 <- numeric(N[ks])
      S1 <- matrix(0, nrow = NP, ncol = N[ks])

      for (i in indices) {
        S2 <- matrix(0, nrow = NP, ncol = NP)

        L <- which(X[ks, 1:N[ks]] >= X[ks, i])

        if(NP == 1){
          ZTL <- ZT[ks,1:NP , L]
          S0[i] <- sum(exp(BZt[L,1])* WGHT[ks, L])
          S1[1, i] <- sum(ZTL*exp(BZt[L,1])* WGHT[ks, L])
          S2 <- S2 + sum(ZTL*ZTL*exp(BZt[L,1])* WGHT[ks, L])

        }else{
          S0[i] <- sum(exp(BZt[L, 1])* WGHT[ks, L])
          S1[1:NP, i] <-  ZT[ks, 1:NP, L] %*% (exp(BZt[L, 1]) *  WGHT[ks, L])
          for(J in 1:NP){
            for(K in 1:NP){
              S2[J, K] <- t(exp(BZt[L, 1])* WGHT[ks, L])%*%(ZT[ks,J , L]*ZT[ks,K , L])
            }
          }
        }
        if (S0[i] > 0.00000001) {

          U <- U + DELTA[ks, i] * (ZT[ks, , i] - S1[, i] / S0[i]) * WGHT[ks, i]


          S1_i <- S1[, i]
          ZT_ks_i <- ZT[ks, , i]
          WGHT_ks_i <- WGHT[ks, i]

          F <- F + DELTA[ks, i] * (S2 / S0[i] - tcrossprod(S1_i) / S0[i]^2) * WGHT_ks_i
          F2 <- F2 + (DELTA[ks, i])^2 * tcrossprod(ZT_ks_i - S1_i / S0[i]) * WGHT_ks_i^2
        }
      }
    }

    BETA <- BETA + solve(F, U)
  }
  invF <- solve(F)
  # C     the variance of \hat\beta(v) is var in the following
  var <- solve(F) %*% F2 %*% solve(F)

  return(list(BETA = BETA, F = invF, var = var))
}

# C     *****************************************************************
# C     *                                                               *
# C     *                    ESTPIPW       function                   *
# C     *                                                               *
# C     * IPW estimate the varying coefficient in the mark PH model     *
# C     *****************************************************************

estpipw <- function(tau, KK, N, TBAND, NP, X, ZT, CENSOR, DELTA, WGHT, beta0) {
  #beta0 is initial value
  BETA = beta0
  var <- matrix(0, nrow = NP, ncol = NP)
  LAMBDA0 <- matrix(NA, nrow = KK, ncol = max(N))

  KL <- 6

  for (KITER in 1:KL) {
    U <- numeric(NP)
    F <- matrix(0, nrow = NP, ncol = NP)
    F2 <- matrix(0, nrow = NP, ncol = NP)

    for (ks in 1:KK) {
      if(NP==1){
        BZt <- as.matrix(ZT[ks, , 1:N[ks]] * BETA, ncol = 1)

      }else{
        BZt <- ZT[ks, , 1:N[ks]] %*% BETA}

      # c*******************************************************************

      mask <- X[ks, 1:N[ks]] <= tau
      indices <- which(mask)
      S0 <- numeric(N[ks])
      S1 <- matrix(0, nrow = NP, ncol = N[ks])

      for (i in indices) {

        S2 <- matrix(0, nrow = NP, ncol = NP)

        L <- which(X[ks, 1:N[ks]] >= X[ks, i])
        if(NP == 1){
          ZTL <- ZT[ks,1:NP , L]
          S0[i] <- sum(exp(BZt[L,1])* WGHT[ks, L])
          S1[1, i] <- sum(ZTL*exp(BZt[L,1])* WGHT[ks, L])
          S2 <- S2 + sum(ZTL*ZTL*exp(BZt[L,1])* WGHT[ks, L])

        }else{
          S0[i] <- sum(exp(BZt[L, 1])* WGHT[ks, L])
          S1[1:NP, i] <-  ZT[ks, 1:NP, L] %*% (exp(BZt[L, 1]) *  WGHT[ks, L])
          for(J in 1:NP){
            for(K in 1:NP){
              S2[J, K] <- t(exp(BZt[L, 1])* WGHT[ks, L])%*%(ZT[ks,J , L]*ZT[ks,K , L])
            }
          }
        }


        if (S0[i] > 0.00000001) {
          U <- U + DELTA[ks, i] * (ZT[ks, , i] - S1[, i] / S0[i]) * WGHT[ks, i]

          F <- F + DELTA[ks, i] * (S2 / S0[i] - tcrossprod(S1[, i]) / S0[i]^2) * WGHT[ks, i]

          F2 <- F2 + DELTA[ks, i]^2 * tcrossprod(ZT[ks, , i] - S1[, i] / S0[i]) * WGHT[ks, i]^2
        }
      }

      # C     CALCULATE THE BASELINE HAZARD FUNCTION AT X(KS,I) AND VSPOT:

      if (KITER == KL) {

        LAMBDA0[ks, 1:N[ks]] <- sapply(1:N[ks], function(i) {
          II <- which((X[ks, 1:N[ks]] <= tau) & (S0[1:N[ks]] > 0.00000001))

          sum(Epanker(X[ks, II], X[ks, i], TBAND, CENSOR[ks, II]) * DELTA[ks, II] * WGHT[ks, II] / S0[II])
        })
      }
    }

    BETA <- BETA + solve(F, U)
  }
  invF <- solve(F)
  # C     the variance of \hat\beta(v) is var in the following
  var <- solve(F) %*% F2 %*% solve(F)

  return(list(BETA = BETA, LAMBDA0 = LAMBDA0, var = var))
}


# C     *****************************************************************
# C     *                                                               *
# C     *                    ESTPAUG       function                   *
# C     *                                                               *
# C     *AIPW estimate the varying coefficient in the mark PH model     *
# C     *****************************************************************

estpaug <- function(tau, tstep, ntgrid, TBAND, KK, N, NP, X, ZT, CENSOR, DELTA, WGHT, DRHOipw, BETA0) {

  mxtgrid <- 100

  var <- matrix(0, nrow = NP, ncol = NP)

  KL <- 6
  BETA <- BETA0

  for (KITER in 1:KL) {
    U <- numeric(NP)
    F <- matrix(0, nrow = NP, ncol = NP)
    F2 <- matrix(0, nrow = NP, ncol = NP)
    S0 <- matrix(0, nrow = KK, ncol = max(N))
    for (ks in 1:KK) {
      if(NP==1){
        BZt <- as.matrix(ZT[ks, , 1:N[ks]] * BETA, ncol = 1)

      }else{
        BZt <- ZT[ks, , 1:N[ks]] %*% BETA
      }

      # c*******************************************************************

      mask <- X[ks, 1:N[ks]] <= tau
      indices <- which(mask)


      S1 <- matrix(0, nrow = NP, ncol = N[ks])

      for (i in indices) {
        S2 <- matrix(0, nrow = NP, ncol = NP)

        L <- which(X[ks, 1:N[ks]] >= X[ks, i])

        if(NP == 1){
          ZTL <- ZT[ks,1:NP , L]
          S0[ks,i] <- sum(exp(BZt[L,1]))
          S1[1, i] <- sum(ZTL*exp(BZt[L,1]))
          S2 <- S2 + sum(ZTL*ZTL*exp(BZt[L,1]))

        }else{
          S0[ks,i] <- sum(exp(BZt[L, 1]))
          S1[1:NP, i] <-  ZT[ks, 1:NP, L] %*% (exp(BZt[L, 1]))
          for(J in 1:NP){
            for(K in 1:NP){
              S2[J, K] <- t(exp(BZt[L, 1]))%*%(ZT[ks,J , L]*ZT[ks,K , L])
            }
          }
        }


        if (S0[ks, i] > 0.00000001) {
          # c  The score function for AIPW estimator:
          U <- U + (ZT[ks, , i] - S1[, i] / S0[ks, i]) *
            (WGHT[ks, i] * DELTA[ks, i] + (1.0 - WGHT[ks, i]) * CENSOR[ks, i] * DRHOipw[ks, i])

          F <- F + (S2 / S0[ks, i] - tcrossprod(S1[, i]) / S0[ks, i]^2) *
            (WGHT[ks, i] * DELTA[ks, i] + (1.0 - WGHT[ks, i]) * CENSOR[ks, i] * DRHOipw[ks, i])

          F2 <- F2 + tcrossprod(ZT[ks, , i] - S1[, i] / S0[ks, i]) *
            (WGHT[ks, i] * DELTA[ks, i] + (1.0 - WGHT[ks, i]) * CENSOR[ks, i] * DRHOipw[ks, i])^2
        }
      }

    }

    BETA <- BETA + solve(F, U)
  }
  invF <- solve(F)
  # C     the variance of \hat\beta(v) is var in the following
  var <- solve(F) %*% F2 %*% solve(F)

  # C     CALCULATE THE BASELINE HAZARD FUNCTION AT X(KS,I) AND VSPOT:
  LAMBDA0 <- sapply(1:KK, function(ks) {
    sapply(1:N[ks], function(i) {
      mask_i <- X[ks, 1:N[ks]] <= tau & (S0[ks, 1:N[ks]] > 0.00000001)
      sum(Epanker(X[ks, mask_i], X[ks, i], TBAND, CENSOR[ks, mask_i]) *
            (WGHT[ks, mask_i] * DELTA[ks, mask_i] / S0[ks, mask_i] +
               (1.0 - WGHT[ks, mask_i]) * CENSOR[ks, mask_i] * DRHOipw[ks, mask_i] / S0[ks, mask_i]))
    })
  })
  LAMBDA0 <- t(LAMBDA0)
  # c     CALCULATE THE BASELINE HAZARD FUNCTION AT X(KS,I) AND VSPOT:
  # c    change to the grid points of (t,v): !!!!!

  LAMBDAk0 <- sapply(1:KK, function(ks) {
    sapply(1:ntgrid, function(Itgrid) {
      tvalue <- tstep * Itgrid

      mask_i <- X[ks, 1:N[ks]] <= tau & (S0[ks, 1:N[ks]] > 0.00000001)

      sum(Epanker(X[ks, mask_i], tvalue, TBAND, CENSOR[ks, mask_i]) *
            (WGHT[ks, mask_i] * DELTA[ks, mask_i] / S0[ks, mask_i] +
               (1.0 - WGHT[ks, mask_i]) * CENSOR[ks, mask_i] * DRHOipw[ks, mask_i] / S0[ks, mask_i]))
    })
  })
  LAMBDAk0 <- t(LAMBDAk0)
  return(list(BETA = BETA, LAMBDA0 = LAMBDA0, LAMBDAk0 = LAMBDAk0, F = invF, var = var))
}


# ******************************************************************
#                       function  GDIST2N
#     This is a modified version of GDIST2 where a different formula
#          is used when doing the change of order of integration
#                   involving the kernel function
#
#          The output of the function GDIST2N is CUMBDIST,
#     which is  the Gaussian multiplier for $n^{1/2}(\hat B_{aug}(v)-B(v))$
#
#          THIS function does not requires smoothing in t.
#                 the results are similar to GDIST1,
#           but GDIST1 may produce nan when TBAND is small
# ******************************************************************

GDIST2N <- function(nvgrid, iskip, zdev, KK, N, NP, X, ZT, betaofv, SigmaInv, S0N, S1N, tempaug, AsigInv) {

  #BU1 <- matrix(0, nrow = NP, ncol = nvgrid)
  CUMBDIST <- matrix(0, nrow = NP, ncol = nvgrid)

  for (ispot in iskip:nvgrid) {
    for (ks in 1:KK) {
      BZt <- matrix(0, nrow = N[ks], ncol = 1)

      if(NP==1){
        BZt <- as.matrix(ZT[ks, , 1:N[ks]] * betaofv[1, ispot], ncol = 1)

      }else{
        BZt <- ZT[ks, , 1:N[ks]] %*% betaofv[, ispot]
      }



      Sx0 <- numeric(N[ks])
      Sx1 <- matrix(0, nrow = NP, ncol = N[ks])

      for (i in 1:N[ks]) {

        indices <- X[ks, 1:N[ks]] >= X[ks, i]
        Sx0[i] <- sum(exp(BZt[indices, 1]) * zdev[ks, indices])

        for (j in 1:NP) {
          Sx1[j, i] <- sum(exp(BZt[indices, 1]) * ZT[ks, j, indices] * zdev[ks, indices])
        }

        if (S0N[ks, i, ispot] > 0.00000001) {
          for (j in 1:NP) {
            for (mspot in iskip:nvgrid) {

              index <- seq((j-1)*NP + 1, NP*j,1)
              indexS1N <- seq((ks-1)*NP+1,ks*NP,1)
              tempBU1 <- sum(AsigInv[index, ispot, mspot] * (ZT[ks, 1:NP, i] - S1N[indexS1N,  i, ispot] / S0N[ks, i, ispot]))
              tempXBU2 <- sum(AsigInv[index , ispot, mspot] * (Sx1[1:NP, i] - Sx0[i] * S1N[indexS1N, i, ispot] / S0N[ks, i, ispot]) / S0N[ks, i, ispot])

              CUMBDIST[j, mspot] <- CUMBDIST[j, mspot] + (tempBU1 * zdev[ks, i] - tempXBU2) * tempaug[ks, i, ispot]
            }
          }
        }
      }
    }
  }

  # nsampa <- sum(N)
  # CUMBDIST <- BU1
  # for (j in 1:NP) {
  #   for (ispot in iskip:nvgrid) {
  #     CUMBDIST[j, ispot] <- BU1[j, ispot]
  #   }
  # }

  return(CUMBDIST)
}


# ******************************************************************
# *  Modified to accommodate stratas
# *
# *                      ESTPLOGIT       function
# *
# ******************************************************************

estplogit <- function(KK, N, NPL, Z, R, D) {

  UL <- numeric(NPL)
  UL2 <- numeric(NPL)
  FL <- matrix(0, nrow = NPL, ncol = NPL)
  FVAR <- array(0, dim = c(KK, NPL, NPL))
  PSI <- matrix(0, nrow = KK, ncol = NPL)
  for (ks in 1:KK) {
    for (iter in 1:20) {
      UL <- numeric(NPL)
      UL2 <- numeric(NPL)
      FL <- matrix(0, nrow = NPL, ncol = NPL)

      TZ <- t(Z[ks, , 1:N[ks]] )%*% (PSI[ks, ])
      P <- exp(TZ) / (1.0 + exp(TZ))

      for (j in 1:NPL) {
        UL[j] <- sum((R[ks, 1:N[ks]] * D[ks, 1:N[ks]] * (1.0 - P) -
                        (1.0 - R[ks, 1:N[ks]]) * D[ks, 1:N[ks]] * P) *
                       Z[ks, j, 1:N[ks]], na.rm = TRUE)
        UL2[j] <- UL[j]

        for (k in 1:NPL) {
          FL[j, k] <- sum((R[ks, 1:N[ks]] * D[ks, 1:N[ks]] +
                             (1.0 - R[ks, 1:N[ks]]) * D[ks, 1:N[ks]]) *
                            P * (1.0 - P) *
                            Z[ks, j, 1:N[ks]] * Z[ks, k, 1:N[ks]], na.rm = TRUE)
        }
      }

      # Solve the system of linear equations using solve()
      #inv <- gaussj(FL, UL)
      PSI[ks, ] <- PSI[ks, ] + solve(FL, UL)
      #PSI[ks, ] <- PSI[ks, ] + inv$B
    }

    FVAR[ks, , ] <- solve(FL)
    #FVAR[ks, , ] <- inv$A
  }

  return(list(PSI = PSI, FVAR = FVAR))
}


#****************************************************************************
# calculate (1/h)*K((t_k-t)/h), K is Epanechnikov  kernel
# K(u)=3/4(1-u^2), -1<u<1
#****************************************************************************

Epanker <- function(tk, tvalue, hband, delt) {
  #tk has the same number of elements with delt as a vector; tvalue is a number
  # Calculate the absolute difference between tk and tvalue
  diff <- abs(tk - tvalue)

  # Create a logical index vector for abs(tk - tvalue) < hband
  idx <- diff < hband
  idx[is.na(idx)] <- FALSE
  # Initialize the result vector
  result <- rep(0,length(tk))

  # Calculate Epanker for elements where abs(tk - tvalue) < hband
  result[idx] <- 3.0 / (4 * hband) * delt[idx] * (1 - (diff[idx] / hband)^2)

  return(result)
}

EpankerV <- function(tk, tvalue, hband, delt) {
  #tk has the same number of elements with delt as a vector; tvalue is a number
  # Calculate the absolute difference between tk and tvalue
  diff <- abs(tk - tvalue)

  # Create a logical index vector for abs(tk - tvalue) < hband
  idx <- diff < hband
  idx[is.na(idx)] <- FALSE
  # Initialize the result vector
  result <- rep(0,length(tk))

  # Calculate Epanker for elements where abs(tk - tvalue) < hband
  result[idx] <- 3.0 / (4 * hband) * delt * (1 - (diff[idx] / hband)^2)

  return(result)
}


EpankerC <- function(tk, tvalue, hband, delt) {
  if (abs(tk - tvalue) < hband) {
    result <- 3.0 / (4 * hband) * delt * (1 - ((tk - tvalue) / hband)^2)
  } else {
    result <- 0.0
  }
  return(result)
}

