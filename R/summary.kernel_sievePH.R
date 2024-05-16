#' Summarizing Nonparametric Kernel-Smoothed Stratified Mark-Specific
#' Proportional Hazards Model Fits
#'
#' \code{summary} method for an object of class \code{kernel_sievePH}.
#'
#' @aliases print.summary.kernel_sievePH
#' @param object an object of class \code{kernel_sievePH}, a result of a call to
#'   \code{\link{kernel_sievePH}} and \code{\link{kernel_sievePHaipw}}.
#' @param contrast a character string specifying the treatment effect parameter
#'   of interest. The default value is \code{"te"} (treatment efficacy); other
#'   options are \code{"hr"} (hazard ratio) and \code{"loghr"} (log hazard
#'   ratio).
#' @param sieveAlternative a character string specifying the alternative
#'   hypothesis for the sieve tests, which can be either \code{"twoSided"}
#'   (default) or \code{"oneSided"}.
#' @param confLevel the confidence level (0.95 by default) of reported
#'   confidence intervals
#'
#' @details \code{print.summary.kernel_sievePH} prints a formatted summary of
#'   results. Inference about coefficients in the kernel-smoothed mark-specific
#'   proportional hazards model is tabulated. Additionally, a summary is
#'   generated
#' from the tests of two relevant null hypotheses: (1) \{\eqn{H_0: HR(v)=1} for
#' all \eqn{v}\}, and (2) \{\eqn{H_0: HR(v)=HR} for all \eqn{v}\}. For the tests
#'   of (2), \code{sieveAlternative} controls the choice of the alternative
#'   hypothesis.
#'
#' @return An object of class \code{summary.kernel_sievePH}, which is a list
#'   with the following components:
#' \itemize{
#' \item \code{estBeta}: a data frame summarizing point estimates and standard
#' errors of the mark-specific coefficients for treatment.
#'
#' \item \code{HRunity.2sided}: a data frame with test statistics (first row) and
#' corresponding p-values (second row) for testing \eqn{H_{10}: HR(v) = 1} vs.
#' \eqn{H_{1a}: HR(v) \neq 1} for any v \eqn{\in [a, b]} (general
#' alternative). \code{TSUP1} is based on an extension of the classic Kolmogorov-Smirnov
#' supremum-based test. \code{Tint1} is a generalization
#' of the integration-based Cramer-von Mises test.
#'
#' \item \code{HRunity.1sided}: a data frame with test statistics (first row) and
#' corresponding p-values (second row) for testing \eqn{H_{10}: HR(v) = 1} vs.
#' \eqn{H_{1m}: HR(v) \leq 1} with strict inequality for some v \eqn{\in [a, b]}
#' (monotone alternative). \code{TSUP1m} is based on an extension of the classic
#' Kolmogorov-Smirnov supremum-based test. \code{Tint1m} is a generalization
#' of the integration-based Cramer-von Mises test.
#'
#'
#' \item \code{HRconstant.2sided}: a data frame with test statistics (first row) and
#' corresponding p-values (second row) for testing \eqn{H_{20}}: HR(v) does not depend on v
#' \eqn{\in [a, b]} vs. \eqn{H_{2a}}: HR depends on v \eqn{\in [a, b]}
#' (general alternative). \code{TSUP2} is based on an extension of the classic
#' Kolmogorov-Smirnov supremum-based test. \code{Tint2} is a generalization
#' of the integration-based Cramer-von Mises test.
#' This component is available if \code{sieveAlternative="twoSided"}.
#'
#' \item \code{HRconstant.1sided}: a data frame with test statistics (first row) and
#' corresponding p-values (second row) for testing \eqn{H_{20}}: HR(v) does not depend on v
#' \eqn{\in [a, b]} vs. \eqn{H_{2m}}: HR increases as v increases \eqn{\in [a, b]}
#' (monotone alternative). \code{TSUP2m} is based on an extension of the classic
#' Kolmogorov-Smirnov supremum-based test. \code{Tint2m} is a generalization
#' of the integration-based Cramer-von Mises test.
#' This component is available if \code{sieveAlternative="oneSided"}.
#'
#' \item \code{te}: a data frame summarizing point and interval estimates of the
#' mark-specific treatment efficacy on the grid of mark values defined by
#' \code{nvgrid} spanning from the minimum and maximum of the mark (available if \code{contrast="te"}).
#' The confidence level is specified by \code{confLevel}.

#' \item \code{hr}: a data frame summarizing point and interval estimates of the
#' mark-specific hazard ratio on the grid of mark values defined by
#' \code{nvgrid} spanning from the minimum and maximum of the mark (available if
#' \code{contrast="hr"}). The confidence level is specified by \code{confLevel}.

#' \item \code{loghr}: a data frame summarizing point and interval estimates of
#' the mark-specific log hazard ratio on the grid of mark values defined by
#' \code{nvgrid} spanning from the minimum and maximum of the mark (available if
#' \code{contrast="loghr"}). The confidence level is specified by
#' \code{confLevel}.


#' }
#'
#' @references
#' Gilbert, P. B. and Sun, Y. (2015). Inferences on relative failure rates in
#' stratified mark-specific proportional hazards models with missing marks, with
#' application to human immunodeficiency virus vaccine efficacy trials.
#' \emph{Journal of the Royal Statistical Society Series C: Applied Statistics},
#' 64(1), 49-73.
#'
#' Sun, Y. and Gilbert, P. B. (2012). Estimation of stratified mark‐specific
#' proportional hazards models with missing marks. \emph{Scandinavian Journal of
#' Statistics}, 39(1), 34-52.
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
#' sfit <- summary(fit)
#' # print the formatted summary
#' sfit
#' # treatment efficacy estimates on the grid
#' sfit$te
#'
#' @seealso \code{\link{kernel_sievePH}}
#'
#' @export
summary.kernel_sievePH <- function(object,
                            contrast = c("te", "hr", "loghr"),
                            sieveAlternative = c("twoSided", "oneSided"), confLevel = 0.95, ...){

  contrast <- match.arg(contrast, choices = c("te", "hr", "loghr"))
  sieveAlternative <- match.arg (sieveAlternative, choices = c("twoSided", "oneSided"))
  # probability level to be used for quantiles in the construction of confidence bounds
  probLevel <- 1 - (1 - confLevel) / 2

  if (!is.null(object$H20)){
    if (sieveAlternative == "oneSided"){
      HRconstant <- object$H20[, c("TSUP2m", "Tint2m")]
    } else {
      HRconstant <- object$H20[, c("TSUP2", "Tint2")]
    }
  }else{
    HRconstant <- NULL
  }

  if(!is.null(object$H10)){
    HRunity.1sided <- object$H10[, c("TSUP1m", "Tint1m")]
    HRunity.2sided <- object$H10[, c("TSUP1", "Tint1")]
  }else{
    HRunity.1sided <- NULL
    HRunity.2sided <- NULL
  }


  estBeta <- object$estBeta

  # linear score lower and upper confidence limits with confidence level specified by 'confLevel'
  betaLB <- estBeta$beta - estBeta$se * qnorm(probLevel)
  betaUB <- estBeta$beta + estBeta$se * qnorm(probLevel)
  betaEst <- estBeta$beta
  # contrast estimates and confidence bounds
  if (contrast == "te") {
    est <- 1 - exp(betaEst)
    lb <- 1 - exp(betaUB)
    ub <- 1 - exp(betaLB)
  } else if (contrast == "hr") {
    est <- exp(betaEst)
    lb <- exp(betaLB)
    ub <- exp(betaUB)
  } else if (contrast == "loghr"){
    est <- betaEst
    lb <- betaLB
    ub <- betaUB
  }

  # assemble the contrast matrix of the output list
  contrastDF <- data.frame(estBeta$mark, est, lb, ub)
  colnames(contrastDF) <- c("mark", switch(contrast, te="TE", hr="HR", loghr="LogHR"), "LB", "UB")

  out <- list(estBeta, HRunity.2sided, HRunity.1sided, HRconstant, contrastDF)
  names(out) <- c("estBeta", "HRunity.2sided", "HRunity.1sided",
                  paste0("HRconstant.", ifelse(sieveAlternative=="oneSided", "1", "2"), "sided"),
                  contrast)

  class(out) <- "summary.kernel_sievePH"
  return(out)
}

#' @rdname summary.kernel_sievePH
#' @param x an object of class \code{summary.kernel_sievePH}, usually a result of a call to \code{summary.kernel_sievePH}
#' @param digits the number of significant digits to use when printing (4 by default)
#' @param ... further arguments passed to or from other methods
#' @export
print.summary.kernel_sievePH <- function(x, digits=4, ...){
  cat("\nCoefficients:\n")
  print(x$estBeta, digits=digits, print.gap=2)
  cat("\n")
  if (!is.null(x$HRunity.2sided) & !is.null(x$HRunity.1sided)) {
    cat("Tests of H0: HR(v) = 1 for all v:\n")
    cat("  Two-sided test:\n")
    cat("    Supremum-based test p-value: ", format(x$HRunity.2sided[2,"TSUP1"], digits=digits, nsmall=digits), "\n", sep="")
    cat("    Integration-based test p-value: ", format(x$HRunity.2sided[2,"Tint1"], digits=digits, nsmall=digits), "\n", sep="")
    cat("  One-sided test:\n")
    cat("    Supremum-based test p-value: ", format(x$HRunity.1sided[2,"TSUP1m"], digits=digits, nsmall=digits), "\n", sep="")
    cat("    Integration-based test p-value: ", format(x$HRunity.1sided[2,"Tint1m"], digits=digits, nsmall=digits), "\n", sep="")
  }

  if (!is.null(x$HRconstant.2sided) | !is.null(x$HRconstant.1sided)) {
    cat ("Tests of H0: HR(v) = HR for all v:\n")
    if (!is.null(x$HRconstant.2sided)) {
      cat("  Supremum-based test p-value: ", format(x$HRconstant.2sided[2,"TSUP2"], digits=digits, nsmall=digits), "\n", sep="")
      cat("  Integration-based test p-value: ", format(x$HRconstant.2sided[2,"Tint2"], digits=digits, nsmall=digits), "\n", sep="")
    } else {
      cat("  Supremum-based test p-value: ", format(x$HRconstant.1sided[2,"TSUP2m"], digits=digits, nsmall=digits), "\n", sep="")
      cat("  Integration-based test p-value: ", format(x$HRconstant.1sided[2,"Tint2m"], digits=digits, nsmall=digits), "\n", sep="")
    }
  }

}
