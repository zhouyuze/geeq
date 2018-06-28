print.geeq <- function(object, ...) {
    Coefs <- matrix(NA,nrow=length(object$beta),ncol=4)
    Coefs[,1] <- c(object$beta)
    Coefs[,2] <- sqrt(diag(object$sandwich))
    Coefs[,3] <- Coefs[,1]/Coefs[,2]
    Coefs[,4] <- round(2*pnorm(abs(Coefs[,3]), lower.tail=F), digits=8)
    colnames(Coefs) <- c("Estimates", "Robust SE", "wald", "p")
    rownames(Coefs) <- c(object$coefnames)
    print(signif(Coefs, digits=4))

    cat("\n Estimated Correlation Parameter: ", signif(object$alpha, digits=4), "\n")
    cat("\n Scale Parameter: ", signif(object$phi, digits=4), "\n")
    cat("\n Correlation Model: ", object$corr.type)
    cat("\n Family Type: ", object$family.type)
    cat("\n Link Type: ", object$link.type, "\n")

    cat("\n Number of GEE iterations:", object$niter)
    cat("\n Number of clusters: ", length(object$cluster.size))
    cat("\n Maximum cluster size: ", max(object$cluster.size), "\n")
}