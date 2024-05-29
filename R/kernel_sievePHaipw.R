#' @useDynLib sievePH, .registration = TRUE
NULL
#' @import np
NULL
#' Nonparametric Kernel-Smoothed Stratified Mark-Specific Proportional Hazards
#' Model with a Univariate Continuous Mark, Missing-at-Random in Some Failures
#'
#' \code{kernel_sievePH} implements estimation methods of Sun and Gilbert (2012)
#' and hypothesis testing methods of Gilbert and Sun (2015) for a mark-specific
#' proportional hazards model accommodating that some failures have a missing
#' mark. The methods allow separate baseline mark-specific hazard functions for
#' different baseline subgroups. Missing marks are handled via augmented IPW
#' (AIPW) approach.
#'
#' @param eventTime a numeric vector specifying the observed right-censored
#'   event time.
#' @param eventInd a numeric vector indicating the event of interest (1 if
#'   event, 0 if right-censored).
#' @param mark a numeric vector specifying a univariate continuous mark subject
#'   to missingness at random. Missing mark values should be set to \code{NA}.
#'   For subjects with \code{eventInd = 0}, the value in \code{mark} should also
#'   be set to \code{NA}.
#' @param tx a numeric vector indicating the treatment group (1 if treatment, 0
#'   if placebo).
#' @param aux a numeric vector specifying a binary or a continuous auxiliary
#'   covariate which may be potentially useful for predicting missingness, i.e,
#'   the probability of missing, and for informing about the distribution of
#'   missing marks. The mark missingness model only requires that the auxiliary
#'   covariates be observed in subjects who experienced the event of interest.
#'   For subjects with \code{eventInd = 0}, the value in \code{aux} may be set
#'   to \code{NA}. If no auxiliary covariate is used, set \code{aux} to the
#'   default of \code{NULL}.
#' @param auxType a character string describing the data type of \code{aux} if
#'   \code{aux} is used. Data types allowed include "binary" and "continuous".
#'   If \code{aux} is not used, \code{auxType} should be set to the default of
#'   \code{NULL}.
#' @param zcov a data frame with one row per subject specifying possibly
#'   time-dependent covariate(s) (not including \code{tx}). If no covariate is
#'   used, \code{zcov} should be set to the default of \code{NULL}.
#' @param strata a numeric vector specifying baseline strata (\code{NULL} by
#'   default). If specified, a separate mark-specific baseline hazard is assumed
#'   for each stratum. It also allows the models of the probability of
#'   complete-case and of the mark distribution to differ across strata.
#' @param formulaPH a one-sided formula object (on the right side of the
#'   \code{~} operator) specifying the linear predictor in the proportional
#'   hazards model. Available variables to be used in the formula include
#'   \code{tx} and variable(s) in \code{zcov}. By default, \code{formulaPH} is
#'   specified as \code{~ tx}.
#' @param formulaMiss a one-sided formula object (on the right side of the
#'   \code{~} operator) specifying the linear predictor in the logistic
#'   regression model used for predicting the probability of observing the mark.
#'   \code{formulaMiss} must be provided for the \code{AIPW} method. Available
#'   variables to be used in the formula include \code{eventTime}, \code{tx},
#'   \code{aux}, and variable(s) in \code{zcov}.
#' @param formulaAux a one-sided formula object (on the right side of the
#'   \code{~} operator) specifying the variables used for estimating the
#'   conditional probability distribution or conditional density of \code{aux}.
#'   If \code{aux} is binary, the formula specifies the linear predictor in a
#'   logistic regression and if \code{aux} is continuous, the formula provides a
#'   symbolic description of variables to be used in kernel conditional density
#'   estimation. \code{formulaAux} is optional for the \code{AIPW} estimation
#'   procedure. Available variables to be used in the formula include
#'   \code{eventTime}, \code{tx}, \code{mark}, and variable(s) in \code{zcov}.
#' @param tau a numeric value specifying the duration of study follow-up period.
#'   Failures beyond \code{tau} are treated as right-censored. There needs to be
#'   at least \eqn{10\%} of subjects (as a rule of thumb) remaining uncensored
#'   by \code{tau} for the estimation to be stable. By default, \code{tau} is
#'   set as the maximum of \code{eventTime}.
#' @param tband a numeric value between 0 and \code{tau} specifying the
#'   bandwidth of the kernel smoothing function over time. By default,
#'   \code{tband} is set as (\code{tau}-min(\code{eventTime}))/5.
#' @param hband a numeric value between 0 and 1 specifying the bandwidth of the
#'   kernel smoothing function over min-max normalized mark which is defined as
#'   \code{(v - min(v)) / (max(v) - min(v))} where \code{v} is the observed
#'   mark. By default, \code{hband} is set as \eqn{4\sigma n^{-1/3}} where
#'   \eqn{\sigma} is the estimated standard deviation of the min-max normalized
#'   marks for uncensored failure times and \eqn{n} is the number of subjects in
#'   the dataset. Larger bandwidths are recommended for higher percentages of
#'   missing marks.
#' @param nvgrid an integer value (100 by default) specifying the number of
#'   equally spaced mark values between the minimum and maximum of the observed
#'   mark for which the treatment effects are evaluated.
#' @param a a numeric value between the minimum and maximum of observed mark
#'   values specifying the lower bound of the range for testing the null
#'   hypotheses \eqn{H_{10}: HR(v) = 1} and \eqn{H_{20}: HR(v)} does not depend
#'   on \eqn{v}, for \eqn{v \in [a, b]}; By default, \code{a} is set as
#'   \code{(max(mark) - min(mark))/nvgrid + min(mark)}.
#' @param b a numeric value between the minimum and maximum of observed mark
#'   specifying the upper bound of the range for testing the null hypotheses
#'   \eqn{H_{10}: HR(v) = 1} and \eqn{H_{20}: HR(v)} does not depend on \eqn{v},
#'   for \eqn{v \in [a, b]}; By default, \code{b} is set as \code{max(mark)}.
#' @param ntgrid an integer value (\code{NULL} by default) specifying the number
#'   of equally spaced time points for which the mark-specific baseline hazard
#'   functions are evaluated. If \code{NULL}, baseline hazard functions are not
#'   evaluated.
#' @param nboot number of bootstrap iterations (500 by default) for simulating
#'   the distributions of the test statistics. If \code{NULL}, the hypotheses
#'   tests are not performed.
#' @param maxit Maximum number of iterations to attempt for convergence in
#'   estimation. The default number is 6.
#' @param seed an integer specifying the random seed for reproducing the test
#'   statistics and p-values. By default, a specific seed is not set.
#'
#' @details \code{kernel_sievePH} analyzes data from a randomized
#'   placebo-controlled trial that evaluates treatment efficacy for a
#'   time-to-event endpoint with a continuous mark. The parameter of interest is
#'   the ratio of the conditional mark-specific hazard functions
#'   (treatment/placebo), which is based on a stratified mark-specific
#'   proportional hazards model. This model assumes no parametric form for the
#'   baseline hazard function nor the treatment effect across different mark
#'   values. For data with missing marks, the estimation procedure leverages
#'   auxiliary predictors of whether the mark is observed and augments the IPW
#'   estimator with auxiliary predictors of the missing mark value.
#'
#' @return An object of class \code{kernel_sievePH} which can be processed by
#'   \code{\link{summary.kernel_sievePH}} to obtain or print a summary of the
#'   results. An object of class \code{kernel_sievePH} is a list containing the
#'   following components:
#' \itemize{
#' \item \code{H10}: a data frame with test statistics (first row) and
#' corresponding p-values (second row) for testing \eqn{H_{10}: HR(v) = 1} for v
#' \eqn{\in [a, b]}. Columns \code{TSUP1} and \code{Tint1} include test
#' statistics and p-values for testing \eqn{H_{10}} vs. \eqn{H_{1a}: HR(v) \neq
#' 1} for any \eqn{v \in [a, b]} (general alternative). Columns \code{TSUP1m}
#' and \code{Tint1m} include test statistics and p-values for testing
#' \eqn{H_{10}} vs. \eqn{H_{1m}: HR(v) \leq 1} with strict inequality for some \eqn{v} in
#' \eqn{[a, b]} (monotone alternative). \code{TSUP1} and \code{TSUP1m} are
#' based on extensions of the classic Kolmogorov-Smirnov supremum-based test.
#' \code{Tint1} and \code{Tint1m} are based on generalizations of the
#' integration-based Cramer-von Mises test. \code{Tint1} and \code{Tint1m}
#' involve integration of deviations over the whole range of the mark. If
#' \code{nboot} is \code{NULL}, \code{H10} is returned as \code{NULL}.
#'
#' \item \code{H20}: a data frame with test statistics (first row) and
#' corresponding p-values (second row) for testing \eqn{H_{20}}: HR(v) does not
#' depend on \eqn{v \in [a, b]}. Columns \code{TSUP2} and \code{Tint2} include
#' test statistics and p-values for testing \eqn{H_{20}} vs. \eqn{H_{2a}}: HR
#' depends on \eqn{v \in [a, b]} (general alternative). Columns \code{TSUP2m}
#' and \code{Tint2m} include test statistics and p-values for testing
#' \eqn{H_{20}} vs. \eqn{H_{2m}}: HR increases as \eqn{v} increases \eqn{\in [a, b]}
#' (monotone alternative). \code{TSUP2} and \code{TSUP2m} are based on
#' extensions of the classic Kolmogorov-Smirnov supremum-based test.
#' \code{Tint2} and \code{Tint2m} are based on generalizations of the
#' integration-based Cramer-von Mises test. \code{Tint2} and \code{Tint2m}
#' involve integration of deviations over the whole range of the mark. If
#' \code{nboot} is \code{NULL}, \code{H20} is returned as \code{NULL}.
#'
#' \item \code{estBeta}: a data frame summarizing point estimates and standard
#' errors of the mark-specific coefficients for treatment at equally-spaced
#' values between the minimum and the maximum of the observed mark values.
#'
#' \item \code{cBproc1}: a data frame containing equally-spaced mark values in
#' the column \code{Mark}, test processes \eqn{Q^{(1)}(v)} for observed data in
#' the column \code{Observed}, and \eqn{Q^{(1)}(v)} for \code{nboot} independent
#' sets of normal samples in the columns S1, S2, \eqn{\cdots}. If
#' \code{nboot} is \code{NULL}, \code{cBproc1} is returned as \code{NULL}.
#'
#' \item \code{cBproc2}: a data frame containing equally-spaced mark values in
#' the column \code{Mark}, test processes \eqn{Q^{(2)}(v)} for observed data in
#' the column \code{Observed}, and \eqn{Q^{(2)}(v)} for \code{nboot} independent
#' sets of normal samples in the columns S1, S2, \eqn{\cdots}. If
#' \code{nboot} is \code{NULL}, \code{cBproc2} is returned as \code{NULL}.
#'
#' \item \code{Lambda0}: an array of dimension K x nvgrid x ntgrid for the
#'   kernel-smoothed baseline hazard function \eqn{\lambda_{0k}, k = 1, \dots,
#'   K} where \eqn{K} is the number of strata. If \code{ntgrid} is \code{NULL}
#' (by default), \code{Lambda0} is returned as \code{NULL}.}
#'
#' @references Gilbert, P. B. and Sun, Y. (2015). Inferences on relative failure
#'   rates in stratified mark-specific proportional hazards models with missing
#'   marks, with application to human immunodeficiency virus vaccine efficacy
#' trials. \emph{Journal of the Royal Statistical Society Series C: Applied
#' Statistics}, 64(1), 49-73.
#'
#'   Sun, Y. and Gilbert, P. B. (2012). Estimation of stratified mark‐specific
#' proportional hazards models with missing marks. \emph{Scandinavian Journal of
#' Statistics}, 39(1), 34-52.
#'
#'   Yang, G., Sun, Y., Qi, L., & Gilbert, P. B. (2017). Estimation of
#'   stratified mark-specific proportional hazards models under two-phase
#'   sampling with application to HIV vaccine efficacy trials. \emph{Statistics
#'   in biosciences}, 9, 259-283.
#'
#' @examples
#' set.seed(20240410)
#' beta <- 2.1
#' gamma <- -1.3
#' n <- 200
#' tx <- rep(0:1, each = n / 2)
#' tm <- c(rexp(n / 2, 0.2), rexp(n / 2, 0.2 * exp(gamma)))
#' cens <- runif(n, 0, 15)
#' eventTime <- pmin(tm, cens, 3)
#' eventInd <- as.numeric(tm <= pmin(cens, 3))
#' alpha <- function(b){ log((1 - exp(-2)) * (b - 2) / (2 * (exp(b - 2) - 1))) }
#' mark0 <- log(1 - (1 - exp(-2)) * runif(n / 2)) / (-2)
#' mark1 <- log(1 + (beta - 2) * (1 - exp(-2)) * runif(n / 2) / (2 * exp(alpha(beta)))) /
#'   (beta - 2)
#' mark <- ifelse(eventInd == 1, c(mark0, mark1), NA)
#' # the true TE(v) curve underlying the data-generating mechanism is:
#' # TE(v) = 1 - exp{alpha(beta) + beta * v + gamma}
#'
#' # a binary auxiliary covariate
#' A <- sapply(exp(-0.5 - 0.2 * mark) / (1 + exp(-0.5 - 0.2 * mark)),
#'             function(p){ ifelse(is.na(p), NA, rbinom(1, 1, p)) })
#' linPred <- 1 + 0.4 * tx - 0.2 * A
#' probs <- exp(linPred) / (1 + exp(linPred))
#' R <- rep(NA, n)
#' while (sum(R, na.rm = TRUE) < 10){
#'   R[eventInd == 1] <- sapply(probs[eventInd == 1],
#'                              function(p){ rbinom(1, 1, p) })
#' }
#' # a missing-at-random mark
#' mark[eventInd == 1] <- ifelse(R[eventInd == 1] == 1, mark[eventInd == 1], NA)
#'
#' # AIPW estimation; auxiliary covariate is used (not required)
#' fit <- kernel_sievePHaipw(eventTime, eventInd, mark, tx, aux = A,
#'                           auxType = "binary", formulaMiss = ~ eventTime,
#'                           formulaAux = ~ eventTime + tx + mark,
#'                           tau = 3, tband = 0.5, hband = 0.3, nvgrid = 20,
#'                           nboot = 20)
#'
#' @importFrom plyr laply
#'
#' @export
kernel_sievePHaipw <- function(eventTime, eventInd, mark, tx, aux = NULL, auxType = NULL, zcov = NULL, strata = NULL, formulaPH = ~ tx, formulaMiss = NULL, formulaAux = NULL,
                               tau = NULL, tband = NULL, hband = NULL, nvgrid = 100, a = NULL, b = NULL, ntgrid = NULL, nboot = 500,  seed = NULL, maxit = 6) {


  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.null(aux)) {
    nauxiliary <- 0
  }else {
    nauxiliary <- 1
  }
  missmark <- 1

  if (min(eventTime) <= 0 | sum(is.na(eventTime)) > 0) {
    stop("eventTime needs to be greater than 0 and have no missing values")
  }

  if (sum(!unique(tx) %in% c(0,1)) > 0) {
    stop("Treatment groups need to take value 0 or 1")
  }

  if (is.null(tau)) {
    tau <- max(eventTime)
  }else {
    tau <- tau
  }
  if (is.null(tband)) {
    tband <- (tau - min(eventTime)) / 5
  } else {
    tband <- tband
  }
  

  nsamp <- NULL
  # kk is the number of strata
  if (is.null(strata)) {
    nsamp[1] <- length(eventTime)
    kk <- 1
  }else {
    stratum <- unique(strata)
    kk <- length(stratum)
    for (i in 1:kk) {
      nsamp[i] <- sum(strata == stratum[i])
    }
  }
  #*****************************************************************************
  ekconst <- 3/5
  #*****************************************************************************
  #     Estimation and Hypothesis tests for Mark PH model with missing marks
  #
  #*****************************************************************************

 # The code inputs mark on its own scale, and re-scales it to 0-1 before doing the analysis.
  vV <- mark
  mn <- min(mark, na.rm = TRUE) #minimum observed mark value
  mx <- max(mark, na.rm = TRUE) #maximum observed mark value
  # Put the marks on the scale 0 to 1:
  vV <- (mark - mn) / (mx - mn)
 
  if (is.null(hband)) {
    hband <- 4*sqrt(var(vV[eventInd == 1], na.rm = TRUE)) * length(eventInd)^{-1/3}
  }else {
    hband <- hband
  }
  
  
  if (is.null(a)) {
    a0 <- 1 / nvgrid
  }else{
    if (a <= mn) {
      stop("a needs to be greater than the minimum of the observed marks")
    }
    a0 <- max(ceiling((a - mn) / (mx - mn) * nvgrid), 1) / nvgrid
  }
  if (is.null(b)) {
    b1 <- 1
  }else{
    if(b > mx) {
      stop("b needs to be no greater than the maximum of the observed marks")
    }
    b1 <- min(floor((b - mn) / (mx - mn) * nvgrid), nvgrid) / nvgrid
  }

  #total sample
  nsampa <- sum(nsamp)
  rnsamp <- nsampa
  #*****************************************************************************
  #                      Set up the smoothing parameters
  #
  #*****************************************************************************

  # Set up the smoothing parameters
  if (!is.null(ntgrid)) {
    tstep <- tau / ntgrid
  }else {
    tstep <- tau / 100
  }

  vstep <- 1 / nvgrid
  iskip <- a0 / vstep

  # Epanechnikov kernel has squared integral = 3/5 !
  ekconst <- 3 / 5

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

  time <- matrix(0, nrow = kk, ncol = max(nsamp))
  txm <- matrix(0, nrow = kk, ncol = max(nsamp))
  censor <- matrix(0, nrow = kk, ncol = max(nsamp))
  markm <- matrix(0, nrow = kk, ncol = max(nsamp))
  Ax <- matrix(0, nrow = kk, ncol = max(nsamp))
  R <- matrix(0, nrow = kk, ncol = max(nsamp))
  if (!is.null(zcov)) {
    dimz <- dim(zcov)[2]
    zcovkk <- array(0, dim = c(kk, dimz, max(nsamp)))
  }
  formulaPHDecomp <- strsplit(strsplit(paste(deparse(formulaPH), collapse = ""), " *[~] *")[[1]], " *[+] *")[[2]]
  formulaPHNew <- as.formula(paste0("eventTime ~ ", paste(formulaPHDecomp, collapse="+")))
  #number of covariates
  ncov <- length(formulaPHDecomp)
  if (!is.null(zcov)) {
    availablePHdf <- data.frame(cbind(tx, zcov))
    colnames(availablePHdf) <- c("tx", colnames(zcov))
  }else {
    availablePHdf <- data.frame("tx" = tx)
    if (ncov >1) {stop("covariates are missing")}
  }
  mfPH <- model.frame(formulaPHNew, availablePHdf)
  #Cox PH regression excludes the intercept
  covartPH <- model.matrix(terms(formulaPHNew), mfPH)[, -1]
  covart <- array(0, dim = c(kk, ncov, max(nsamp)))
  if (kk == 1) {
    n_ks <- nsamp[1]
    ks <- 1
    if (ncov == 1) {
      covart[ks, 1:ncov, 1:n_ks] <- covartPH[1:n_ks]
    }else {
      covart[ks, 1:ncov, 1:n_ks] <- t(covartPH[1:n_ks, 1:ncov])
    }
    if (!is.null(zcov)) {
      zcovkk[ks, 1:dimz, 1:n_ks] <- t(zcov[1:n_ks,1:dimz])
    }
    time[ks, 1:n_ks] <- eventTime[1:n_ks]
    censor[ks, 1:n_ks] <- eventInd[1:n_ks]
    txm[ks, 1:n_ks] <- tx[1:n_ks]
    if (!is.null(aux)) {
      Ax[ks, 1:n_ks] <- aux[1:n_ks]
    }
    R[ks, 1:n_ks] <- 1 - is.na(mark[1:n_ks])
    # Important For all eventInd ==0, R is set as 1
    R[ks,eventInd == 0] <- 1
    markm[ks, 1:n_ks] <- vV[1:n_ks]
  }else{
    for (ks in 1:kk) {
      n_ks <- nsamp[ks]
      indice_ks <- (strata == stratum[ks])
      if (ncov == 1) {
        covart[ks, 1:ncov, 1:n_ks] <- covartPH[indice_ks]
      }else {
        covart[ks, 1:ncov, 1:n_ks] <- t(covartPH[indice_ks, 1:ncov])
      }
      if (!is.null(zcov)) {
        zcovkk[ks, 1:dimz, 1:n_ks] <- t(zcov[indice_ks, 1:dimz])
      }
      time[ks, 1:n_ks] <- eventTime[indice_ks]
      censor[ks, 1:n_ks] <- eventInd[indice_ks]
      markm[ks, 1:n_ks] <- vV[indice_ks]
      txm[ks, 1:n_ks] <- tx[indice_ks]
      if (!is.null(aux)) {
        Ax[ks, 1:n_ks] <- aux[indice_ks]
      }
      R[ks, 1:n_ks] <- 1 - is.na(mark[indice_ks])
      # Important For all eventInd ==0, R is set as 1
      R[ks, censor[ks, 1:n_ks] == 0] <- 1

    }
  }

  # Calculate sumdelta and Rsum
  sumdelta <- sum(censor)
  Rsum <- sum(R * censor)

  # Calculate the probability of nonmissing marks.
  # wght(ks,i) can be calculated using an external program
  # setup the covariates for logit regression for R:

  n_max <- max(nsamp)
  # Handle missmark
  if (missmark == 0) {
    R <- matrix(1, nrow = kk, ncol = n_max)
    wght <- R
  } else {
    if (is.null(formulaMiss)) {
      stop("`formulaMiss` must be provided")}
    formulaMissDecomp <- strsplit(strsplit(paste(deparse(formulaMiss), collapse = ""), " *[~] *")[[1]], " *[+] *")[[2]]
    formulaMissNew <- as.formula(paste0("R ~ ", paste(formulaMissDecomp, collapse="+")))
    # Define ncovR
    ncovR <- length(formulaMissDecomp) + 1
    # Create covartR array
    covartR <- array(0, dim = c(kk, ncovR, n_max))
    for (ks in 1:kk) {
      if (!is.null(aux)) {
        availabledf <- data.frame("R" = R[ks, 1:nsamp[ks]], "eventTime" = time[ks, 1:nsamp[ks]],
        "aux" = Ax[ks, 1:nsamp[ks]], "tx" = txm[ks, 1:nsamp[ks]])
        if (!is.null(zcov)) {
          availabledf <- cbind(availabledf, t(zcovkk[ks, 1:dimz, 1:nsamp[ks]]))
          colnames(availabledf) <- c("R", "eventTime", "aux", "tx", colnames(zcov))
        }else {
          colnames(availabledf) <- c("R", "eventTime", "aux", "tx")
        }
        #this is to avoid mf below remove all rows with NA. The logistic model does not use non-cases so it doesn't care the values of aux for them
        availabledf$aux[is.na(availabledf$aux)] <- 0
      }else{
        availabledf <- data.frame("R" = R[ks, 1:nsamp[ks]], "eventTime" = time[ks, 1:nsamp[ks]], "tx" = txm[ks, 1:nsamp[ks]])
        if (!is.null(zcov)) {
          availabledf <- cbind(availabledf, t(zcovkk[ks, 1:dimz, 1:nsamp[ks]]))
          colnames(availabledf) <- c("R", "eventTime", "tx", colnames(zcov))
        }else {
          colnames(availabledf) <- c("R", "eventTime", "tx")
        }
        if ("aux" %in% formulaMissDecomp) {
          stop("`aux` is not available")
        }
      }

      mf <- model.frame(formulaMissNew, availabledf)
      covartR[ks,1:ncovR,1:nsamp[ks]] <- t(model.matrix(terms(formulaMissNew ), mf))

    }

    psiRse <- matrix(0, nrow = kk, ncol = ncovR)

    estplogitans <- estplogit(kk, nsamp, ncovR, covartR, R, censor)
    psiR <- estplogitans$PSI
    psiRvar <- estplogitans$FVAR

    for (ks in 1:kk) {
      psiRse[ks, ] <- sqrt(diag(psiRvar[ks, , ]))
    }


    #logistic regression
    wght <- matrix(0, nrow = kk, ncol = n_max)
    for (ks in 1:kk) {
     ei <- covartR[ks, 1, ] * psiR[ks, 1] + covartR[ks, 2, ] * psiR[ks, 2]
     phi <- exp(ei) / (1 + exp(ei))
     wght[ks, ] <- R[ks, ] / (censor[ks, ] * phi + 1 - censor[ks, ])
   }
  }


  # the end of calculating wght(ks,i)

  # Calculate the pdf G(ks,i,ispot) for the second stage estimator, aipw
  # IF use (Nauxiliary.eq.1),
  # the pdf G(ks,i,ispot) can also be calculated using an external program
  # The relevant code is surrounded by the comments AUXILIARY CODE

  # AUXILIARY CODE
  # setup the covariates calculate g(a|w):

  if (!is.null(aux)) {
    if (is.null(formulaAux)) {
      stop("`formulaAux` must be provided")
      }
    formulaAuxDecomp <- strsplit(strsplit(paste(deparse(formulaAux), collapse = ""), " *[~] *")[[1]], " *[+] *")[[2]]
    formulaAuxNew <- as.formula(paste0("aux ~ ", paste(formulaAuxDecomp, collapse="+")))
    if (auxType == "binary") {
      ncovG <- length(formulaAuxDecomp) + 1
    }else{
      ncovG <- length(formulaAuxDecomp)
    }
    covartG <- array(0, dim = c(kk, ncovG, n_max))
    for (ks in 1:kk) {
      availabledf <- data.frame("aux" = Ax[ks, 1:nsamp[ks]], "eventTime" = time[ks, 1:nsamp[ks]],
                                "tx" = txm[ks, 1:nsamp[ks]], "mark" = markm[ks, 1:nsamp[ks]])
      #this is to avoid mf below remove all rows with NA (in aux). covartG does not include aux
      #rows with missing mark would not enter estimation
      availabledf$aux[is.na(availabledf$aux)] <- 0
      availabledf$mark[is.na(availabledf$mark)] <- 0
      if (!is.null(zcov)) {
        availabledf <- cbind(availabledf, t(zcovkk[ks, 1:dimz, 1:nsamp[ks]]))
        colnames(availabledf) <- c("aux", "eventTime", "aux", "tx", "mark", colnames(zcov))
      }else{
        colnames(availabledf) <- c("aux", "eventTime", "tx", "mark")
      }
      mf <- model.frame(formulaAuxNew, availabledf)
      if (auxType == "binary") {
        covartG[ks, 1:ncovG, 1:nsamp[ks]] <- t(model.matrix(terms(formulaAuxNew), mf))
      }else {
        #intercept is removed for kernel conditional density estimation
        covartG[ks, 1:ncovG, 1:nsamp[ks]] <- t(model.matrix(terms(formulaAuxNew), mf)[,-1])
        covartGvarNames <- colnames(mf)[-1]
      }

    }

  }
  # Compute cenG; the logistic regression is only based on cases with observed marks
  cenG <- censor * R

  if (nauxiliary == 1) {
    G <- array(1, dim = c(kk, n_max, nvgrid))
    if(auxType == "binary") {
      estplogitGans <- estplogit(kk, nsamp, ncovG, covartG, Ax, cenG)
      psiG <- estplogitGans$PSI
      for (ks in 1:kk) {
        for (i in 1:nsamp[ks]) {
          ei <- 0
          for (j in 1:(ncovG - 1)) {
            ei <- ei + covartG[ks, j, i] * psiG[ks, j]
          }
          vspot <- (1:nvgrid) * vstep
          phi <- exp(ei + psiG[ks, ncovG] * vspot) / (1 + exp(ei + psiG[ks, ncovG] * vspot))
          G [ks, i, 1:nvgrid] <- phi * Ax[ks, i] + (1 - phi) * (1 - Ax[ks, i])
        }
      }
    }else{
      #kernel conditional density estimation
      for (ks in 1:kk) {
        data <- data.frame(cbind(Ax[ks, 1:nsamp[ks]], t(covartG[ks, 1:ncovG,1:nsamp[ks]])))
        colnames(data) <- c("aux", covartGvarNames)
        #use cases with observed mark to conduct the estimation
        cenGind <- which(cenG[ks,1:nsamp[ks]] == 1)
        subdata <- data[cenGind, ]
        xdat <- subdata[, -1]
        ydat <- subdata[, 1]
        bw <- npcdensbw(xdat = xdat,ydat = ydat, cxkertype="epanechnikov", cykertype="epanechnikov")
        vspot <- (1:nvgrid)*vstep
        fhat <- npcdens(bws = bw, txdat = xdat, tydat = ydat)
        indices <- which(!is.na(Ax[ks, 1:nsamp[ks]]))
        for (ispot in 1:nvgrid) {
          newdata <- data
          newdata$mark <- vspot[ispot]
          newdata <- newdata[indices, ]
          exdat <- newdata[, -1]
          eydat <- newdata[, 1]
          G[ks, indices, ispot] <- predict(fhat, exdat = exdat, eydat = eydat)
          #the original code would have G be NA for aux that is NA
          #in code below, only G for noncases are used
        }
      }
    }
  } else {
    G <- array(1, dim = c(kk, n_max, nvgrid))
  }

  # the end of calculating G(ks,i,ispot).
  # AUXILIARY CODE

  coxfit <- estpvry(tau, kk, nsamp, ncov, time, covart, censor)
  betacox <- coxfit$BETA
  varcox <- coxfit$var
  betacom <- matrix(0, ncol = nvgrid, nrow = ncov)
  secom <- matrix(0, ncol = nvgrid, nrow = ncov)
  betaipw <- matrix(0, ncol = nvgrid, nrow = ncov)
  seipw <- matrix(0, ncol = nvgrid, nrow = ncov)
  LAMBDA0ipw <- array(0, dim = c(kk, max(nsamp), nvgrid))
  betaaug <- matrix(0, ncol = nvgrid, nrow = ncov)
  seaug <- matrix(0, ncol = nvgrid, nrow = ncov)
  LAMBDA0aug <- array(0, dim = c(kk, max(nsamp), nvgrid))
  betaofv <- array(0, dim = c(ncov, nvgrid))
  SigmaInv <- array(0, dim = c(ncov, ncov, nvgrid))

  # baseline hazard functions:
  Lamktv0 <- array(0, dim = c(kk, ifelse(is.null(ntgrid), 1, ntgrid), nvgrid))
  if (!is.null(ntgrid)) {
    estBaseLamInd <- 1
  }else{
    estBaseLamInd <- 0
  }
  ###################################################################
  for (ispot in 1:nvgrid) {
    vspot <- ispot * vstep
    deltak <- matrix(0, nrow = kk, ncol = max(nsamp))
    for (ks in 1:kk) {
      deltak[ks, ] <- Epanker(markm[ks, ], vspot, hband, censor[ks, ])
    }
    # Complete data likelihood estimator: complete
    # Set the weight wght(ks, i) = R(ks, i) to get the complete data likelihood estimator
    comfit <- estpcomcplusplus(tau, kk, nsamp, ncov, time, covart, deltak, R, betacox)
    betacom[, ispot] <- comfit$BETA
    secom[, ispot] <- sqrt(diag(comfit$var))


    # The inverse probability weighted estimator: ipw
    ipwfit <- estpipwcplusplus(tau, tstep, ifelse(is.null(ntgrid), 1, ntgrid), tband, kk, nsamp,
                                ncov, time, covart, censor, deltak, wght, betacom[,ispot], maxit, 0)

    betaipw[, ispot] <- ipwfit$BETA
    seipw[, ispot] <- sqrt(diag(ipwfit$var))
    LAMBDA0ipw[, , ispot] <- ipwfit$LAMBDA0
  }



  #   ********************************************************************
  #                         second stage estimator:   aipw
  #   ********************************************************************

    # Second stage estimator: aipw
    CLAMBDAIPW <- array(0, dim = c(kk, max(nsamp)))
    LAMBDAIPW <- array(0, dim = c(kk, max(nsamp), nvgrid))

    LAMBDAUG <- array(0, dim = c(kk, max(nsamp), nvgrid))
    CLAMBDAUG <- array(0, dim = c(kk, max(nsamp)))

    for (ispot in 1:nvgrid) {
      vspot <- ispot * vstep
      # CALCULATE CLAMBDA(KS,I)= $\int_0^1 \lambda_k(t,u|z) h_k(a|t,u,z) du$ AT $W_{ki}=(T_{ki},Z_{ki})$:
      for (ks in 1:kk) {
        valid_indices <- which((time[ks, 1:nsamp[ks]] <= tau) & (censor[ks,1:nsamp[ks] ] > 0.5))

        #vectorize over loop I = 1, NSAMP(ks)
        if (any(valid_indices)) {
          BZIPW <- array(0, dim = c(max(nsamp),1))
          if(ncov == 1){
            BZIPW[valid_indices, 1] <- covart[ks, 1, valid_indices] * betaipw[1, ispot]
          }else{
            BZIPW[valid_indices, 1] <- t(covart[ks, , valid_indices]) %*% betaipw[, ispot] #summing over J
          }

          #the original fortran code had LAMBDA0 for nvgrid only
          LAMBDAIPW[ks, valid_indices, ispot] <- LAMBDA0ipw[ks, valid_indices,ispot] * exp(BZIPW[valid_indices,1])
          CLAMBDAIPW[ks, valid_indices] <- CLAMBDAIPW[ks, valid_indices] + LAMBDA0ipw[ks, valid_indices,ispot] * exp(BZIPW[valid_indices,1]) * G[ks,valid_indices,ispot]*vstep
        }
      }
    }
# Note: LAMBDAIPW is lambda_k(t,u|z) in equation 11, Sun and Gilbert 2012
#      CLAMBDAIPW is the int_0^1 lambda_k(t,u|z)g_k(a|t,u,z)du in equation 11
########################################################
# CALCULATE CKLAMBDA(KS,I, u)= $ K_h(u-v) \lambda_k(t,u|z) h_k(a|t,u,z) du$ AT $W_{ki}=(T_{ki},Z_{ki})$:
# K_h(u-v)du is
#DRHOipw is the derivative of rho_k^{ipw}(v,w_i) w.r.t v in equation 10.
    CKLAMBDAIPW <- array(0, dim = c(kk, max(nsamp), nvgrid))
    for (ks in 1:kk) {
      for (i in 1:nsamp[ks]) {
        vspot <- (1:nvgrid) * vstep
        kernel_matrix <- laply(vspot, function(x) {
          z <- EpankerV(vspot, x, hband, censor[ks, i])
        })
        CKLAMBDAIPW[ks, i, ] <- kernel_matrix %*% (G[ks, i, ]*LAMBDAIPW[ks, i, ]) * vstep
      }
    }

    CLAMBDAUG <- matrix(0, nrow = kk, ncol = max(nsamp))

    CUMB1 <- rep(0, nvgrid)

    vspot <- (1:nvgrid) * vstep

    for (ispot in 1:nvgrid) {
      vspot <- ispot * vstep
      DRHOipw <- matrix(0, ncol = max(nsamp), nrow = kk)
      deltak <- matrix(0, ncol = max(nsamp), nrow = kk)
      for (ks in 1:kk) {
        valid_indices <- which((time[ks, 1:nsamp[ks]] <= tau) & (censor[ks, 1:nsamp[ks]] > 0.5) & (CLAMBDAIPW[ks, 1:nsamp[ks]] > 0.00000001))
        deltak[ks, ] <- Epanker(markm[ks, ], vspot, hband, censor[ks, ])
        DRHOipw[ks, valid_indices] <- CKLAMBDAIPW[ks, valid_indices, ispot] / CLAMBDAIPW[ks, valid_indices]
        DRHOipw[ks, !valid_indices] <- 0.0
      }

      # initial value
      betaaug0 <- betaipw[, ispot]
      estpaug_result = estpaugcplusplus(tau, tstep, ifelse(is.null(ntgrid), 1, ntgrid), tband, kk, nsamp, ncov, time, covart,
                                        censor, deltak, wght, DRHOipw, betaaug0, maxit,
                                        ifelse(estBaseLamInd ==1, 1, 0))
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
      LAMBDA0aug <- estpaug_result$LAMBDA0
      # CALCULATE CLAMBDA(KS,I)= $\int_0^1 \lambda_k(t,u|z) h_k(a|t,u,z) du$ AT $W_{ki}=(T_{ki},Z_{ki})$:
      for (ks in 1:kk) {
        valid_indices <- which((time[ks, 1:nsamp[ks]] <= tau) & (censor[ks, 1:nsamp[ks]] > 0.5))
        if(ncov == 1){
          BZaug <- as.matrix(covart[ks, , valid_indices] * betaaug[1, ispot],ncol=1)
        }else{
          BZaug <- t(covart[ks, , valid_indices]) %*% betaaug[, ispot]
        }
        LAMBDAUG[ks, valid_indices, ispot] <- LAMBDA0aug[ks, valid_indices] * exp(BZaug)
        CLAMBDAUG[ks, valid_indices] <- CLAMBDAUG[ks, valid_indices] + G[ks, valid_indices, ispot] * LAMBDA0aug[ks, valid_indices] * exp(BZaug) * vstep
      }
    }


    if(!is.null(nboot)){
      # **************************************************************************
      # calculate the quantities needed for GDIST2N
      # **************************************************************************
      AsigInv <- array(0, dim = c(ncov * ncov, nvgrid, nvgrid))
      tempaug <- array(0, dim = c(kk, max(nsamp), nvgrid))
      for (ispot in iskip:nvgrid) {
        va <- vstep * (ispot - 1.0)
        vb <- vstep * ispot
        vspot <- ispot * vstep
        sigtemp <- matrix(0.0, nrow = ncov, ncol = ncov)
        for(mspot in iskip:nvgrid){
          vmspot <- mspot * vstep
          sigtemp <- sigtemp + SigmaInv[, , mspot] * Epanker(vspot, vmspot, hband, 1.0) * vstep
          AsigInv[, ispot, mspot] <- as.vector(sigtemp[1:ncov, 1:ncov])
        }
        DRHOaug <- matrix(0, nrow = kk, ncol = max(nsamp))
        for(ks in 1:kk){
          valid_indices <- which((time[ks,] <= tau) & (censor[ks,] > 0.5) & (CLAMBDAUG[ks,] > 0.00000001))
          DRHOaug[ks, valid_indices] <- G[ks, valid_indices, ispot] * LAMBDAUG[ks, valid_indices, ispot] / CLAMBDAUG[ks, valid_indices]
          DRHOaug[ks, !valid_indices] <- 0.0
          temp1 <- (1.0 - wght[ks, ]) * DRHOaug[ks, ] * vstep
          temp1[!(time[ks, ] <= tau) | !(censor[ks, ] > 0.5)] <- 0.0
          temp2 <- wght[ks, ]
          temp2[!(time[ks, ] <= tau) | !(markm[ks, ] > va) | !(markm[ks, ] <= vb) | !(censor[ks, ] > 0.5)] <- 0.0
          tempaug[ks, 1:nsamp[kk], ispot] <- temp1 + temp2
        }
      }

      #############################################
      S0N <- array(0, dim = c(kk, max(nsamp),nvgrid))
      S1N <- array(0, dim = c(kk * ncov, max(nsamp), nvgrid))

      for (ispot in iskip:nvgrid) {
        for (ks in 1:kk) {
          # for time independent covariates
          BZt <- matrix(0, nrow = nsamp[ks], ncol = 1)
          if(ncov == 1){
            BZt[1:nsamp[ks], 1] <- as.matrix(covart[ks, , 1:nsamp[ks]] * betaofv[1, ispot], ncol = 1)
          }else{
            BZt <- t(covart[ks, , 1:nsamp[ks]]) %*% betaofv[, ispot]
          }
          S0 <- rep(0, nsamp[ks])
          S1 <- matrix(0, nrow = ncov, ncol = nsamp[ks])
          for (i in 1:nsamp[ks]) {
            valid_indices <- which(time[ks, 1:nsamp[ks]] >= time[ks, i])
            S0[i] <- sum(exp(BZt[valid_indices, 1]))
            for(J in 1:ncov){
              S1[J,i] <- sum(covart[ks,J , valid_indices]*exp(BZt[valid_indices, 1]))
            }
          }
          S0N[ks, 1:nsamp[ks], ispot] <- S0
          arrayseq <- seq((ks-1)*ncov +1, (ks*ncov), 1)
          # browser()
          for (a in arrayseq) {
            S1N[a, 1:nsamp[ks], ispot] <- S1[a, ]
          }
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
      BootDist <- matrix(0, nrow = nboot, ncol = nvgrid)
      for (iboot in 1:nboot) {
        zdev <- matrix(rnorm(kk * max(nsamp)), nrow = kk, ncol = max(nsamp))
        #CUMBDIST is simulated B(v)
        CUMBDIST <- GDIST2Ncplusplus(nvgrid, iskip, zdev, kk, nsamp, ncov, time, covart,
                                     betaofv, SigmaInv, S0N, S1N, tempaug, AsigInv)
        CUMB1x[iskip:nvgrid] <- CUMB1x[iskip:nvgrid] + CUMBDIST[1, iskip:nvgrid]
        CUMB1se[iskip:nvgrid] <- CUMB1se[iskip:nvgrid] + CUMBDIST[1, iskip:nvgrid]^2.0
        BootDist[iboot, ] <- CUMBDIST[1, ]
      }
      CUMB1x[iskip:nvgrid] <- CUMB1x[iskip:nvgrid] / nboot
      CUMB1se[iskip:nvgrid] <- sqrt((CUMB1se[iskip:nvgrid] - (CUMB1x[iskip:nvgrid])^2.0 * nboot) / nboot)

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
      markgrid <- (1:nvgrid) * vstep
      markgrid_original_scale <- markgrid*(mx-mn)+mn
      cBproc1.df <- data.frame(cbind(markgrid_original_scale[(iskipa1 + 1):nvgrid], vstep*((iskipa1 + 1):nvgrid),
                                     cBproc1[(iskipa1 + 1):nvgrid], cBproc1s[(iskipa1 + 1):nvgrid,]))
      colnames(cBproc1.df) <- c("Mark", "Standardized mark","Observed", paste0("S", 1:nboot))
      cBproc2.df <- data.frame(cbind(markgrid_original_scale[(iskipa2):nvgrid], vstep*((iskipa2):nvgrid),
                                     cBproc2[iskipa2:nvgrid], cBproc2s[iskipa2:nvgrid, ]))
      colnames(cBproc2.df) <- c("Mark", "Standardized mark", "Observed", paste0("S", 1:nboot))
      cBproc1 <- cBproc1.df[, -2]
      cBproc2 <- cBproc2.df[, -2]
      H10 <- ans10
      H20 <- ans20

    }else{
      cBproc1 <- NULL
      cBproc2 <- NULL
      H10 <- NULL
      H20 <- NULL
      markgrid <- (1:nvgrid) * vstep
      markgrid_original_scale <- markgrid*(mx-mn)+mn
    }

    estBeta <- data.frame("mark" = markgrid_original_scale,
                          "betaaug" = betaaug[1, ],
                          "seaug" = seaug[1, ])
    colnames(estBeta) <- c("mark", "beta", "se")

    if(estBaseLamInd == 1){
      Lambda0 <- Lamktv0
    }else{
      Lambda0 <- NULL
    }

    out <- list("cBproc1" = cBproc1,
                "cBproc2" = cBproc2,
                "H10" = H10,
                "H20" = H20,
                "estBeta" = estBeta,
                "Lambda0" = Lambda0)

  class(out) <- "kernel_sievePH"
  return(out)

}




# C     *****************************************************************
# C
# C                         _ESTPVRY       function_
# C
# C     _estimate the time-varying coefficient in the cox model_
# C     *****************************************************************

estpvry <- function(tau, KK, N, NP, X, ZT, DELTA) {

  BETA <- rep(0, NP)
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
        BZt <- t(ZT[ks, , 1:N[ks]]) %*% BETA
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
          ZTL <- ZT[ks, 1:NP , L]
          S0[i] <- sum(exp(BZt[L,1]))
          S1[1, i] <- sum(ZTL * exp(BZt[L, 1]))
          S2 <- S2 + sum(ZTL * ZTL*exp(BZt[L, 1]))
        }else{
          S0[i] <- sum(exp(BZt[L, 1]))
          S1[1:NP, i] <-  ZT[ks, 1:NP, L] %*% exp(BZt[L, 1])
          for(J in 1:NP){
            for(K in 1:NP){
              S2[J, K] <- t(exp(BZt[L, 1])) %*% (ZT[ks, J , L] * ZT[ks, K , L])
            }
          }
        }
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
    for (iter in 1:6) {
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
      PSI[ks, ] <- PSI[ks, ] + solve(FL, UL)
    }
    FVAR[ks, , ] <- solve(FL)
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


