print.summary.sievePHaipw <- function(x, digits=4, ...){
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
