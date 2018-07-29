print.geeq <- function(object, ...) {
    Coefs <- matrix(NA,nrow=length(object$beta),ncol=4)
    Coefs[,1] <- c(object$beta)
    Coefs[,2] <- sqrt(diag(object$variance))
    Coefs[,3] <- Coefs[,1]/Coefs[,2]
    Coefs[,4] <- round(2*pnorm(abs(Coefs[,3]), lower.tail=F), digits=8)
    colnames(Coefs) <- c("Estimates", "Model SE", "wald", "p")
    rownames(Coefs) <- names(object$beta)
    print(signif(Coefs, digits=4))

    cat("\n Correlation Model: ", object$corr.type)
    if (object$method == "gee") {
      cat("\n Correlation Matrix:\n")
      print(signif(object$correlation, digits=4))
    }
    cat("\n Scale Parameter: ", signif(object$phi, digits=4), "\n")

    cat("\n Null deviance: ", signif(object$null.deviance, digits=1), "on ",
      object$df.null, " degrees of freedom")
    cat("\n Residual deviance: ", signif(object$deviance, digits=1), "on ",
      object$df.residual, " degrees of freedom\n")

    if (object$method == "gee") {
      cat(" QIC: ", object$QIC, "\n")
    } else if (object$method == "qif") {
      cat(" Statistics: \n")
      print(signif(object$statistics, digits=4))
    }

    cat("\n Converged: ", object$converged)
    cat("\n Number of iterations: ", object$niter)
    cat("\n Number of clusters: ", length(object$cluster.size))
    cat("\n Maximum cluster size: ", max(object$cluster.size), "\n")
}
