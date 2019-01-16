print.summary.sievePH <- function(x,...){
  cat("\nCoefficients:\n")
  print(coef, digits=4, print.gap=2)
  cat("\n")
  cat("Two-sided likelihood ratio test of H0: HR(v)=1 p-value: ", x$pLR.HRunity.2sided, "\n", sep="")
  cat("Two-sided Wald test of H0: HR(v)=1 p-value: ", x$pWald.HRunity.2sided, "\n", sep="")
  cat("One-sided weighted Wald-type test of H0: HR(v)=1 p-value: ", x$pWtWald.HRunity.1sided, "\n", sep="")
  if (is.null(x$pLR.HRconstant.2sided)) {
    cat("One-sided likelihood ratio test of H0: HR(v)=HR p-value: ", x$pLR.HRconstant.1sided, "\n", sep="")
    cat("One-sided Wald test of H0: HR(v)=HR p-value: ", x$pWald.HRconstant.1sided, "\n", sep="")
  } else {
    cat("Two-sided likelihood ratio test of H0: HR(v)=HR p-value: ", x$pLR.HRconstant.2sided, "\n", sep="")
    cat("Two-sided Wald test of H0: HR(v)=HR p-value: ", x$pWald.HRconstant.2sided, "\n", sep="")
  }
}