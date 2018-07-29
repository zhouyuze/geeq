library(bindata)
library(geeM)
generatedata <- function(n, t) {
  x1 <- rbinom(n * t, 1, 0.5)
  x2 <- rep(-1:(t-2), n)
  x3 <- runif(n * t, -1, 1)
  x4 <- runif(n * t, -1, 1)

  beta <- c(-0.25, -0.25)
  eta <- 0.25 + cbind(x1, x2) %*% beta
  family <- binomial()
  mu <- family$linkinv(eta)

  rho <- 0.5
  cor <- matrix(rho, t, t)
  diag(cor) <- 1

  y <- vector("numeric", t*n)
  for (i in 1:n) {
    y[((i-1)*t+1) : (i * t)] <- rmvbin(1, mu[((i-1)*t+1) : (i * t)], bincorr=cor)
  }
  id <- rep(1:n, each=t)
  data.frame(y=y, x1=x1, x2=x2, x3=x3, x4=x4, id=id)
}

QIC <- function(y, X, family, beta, sandwich, scale = 1, wt = NULL) {
  if (is.null(wt)) {
    wt = rep(1, length(y))
  }
  family <- family()
  eta <- X %*% beta
  mu <- family$linkinv(eta)
  deriv <- family$mu.eta(eta)
  var <- family$variance(mu)

  p <- dim(X)[2]
  Omega <- matrix(0, ncol = p, nrow = p)
  for (i in 1:p) {
    for (j in 1:p) {
      Omega[i, j] <- sum(X[, i] * X[, j] * deriv * deriv  / var) / scale
    }
  }
  dev.resids <- family$dev.resids
  likelyhood <- sum(dev.resids(y, mu, wt)) / scale
  penalty <- 2 * sum(diag((Omega %*% sandwich)))
  return(list(ll=likelyhood, penalty=penalty, QIC=likelyhood+penalty))
}

# n is the number of clusters, each cluster has 3 observations.
variable_select <- function(n, repeat_time, cor_str) {
  qic <- rep(0, 5)
  select <- matrix(0, ncol=5, nrow=1)
  for (i in 1:repeat_time) {
    data <- generatedata(n, 3)
    r <- geem(y~x1, id, data=data, family=binomial, corstr=cor_str)
    qic[1] <- QIC(r$y, r$X, binomial, r$beta, r$var, r$phi)$QIC
    r <- geem(y~x1+x2, id, data=data, family=binomial, corstr=cor_str)
    qic[2] <- QIC(r$y, r$X, binomial, r$beta, r$var, r$phi)$QIC
    r <- geem(y~x1+x3, id, data=data, family=binomial, corstr=cor_str)
    qic[3] <- QIC(r$y, r$X, binomial, r$beta, r$var, r$phi)$QIC
    r <- geem(y~x1+x2+x3, id, data=data, family=binomial, corstr=cor_str)
    qic[4] <- QIC(r$y, r$X, binomial, r$beta, r$var, r$phi)$QIC
    r <- geem(y~x1+x2+x3+x4, id, data=data, family=binomial, corstr=cor_str)
    qic[5] <- QIC(r$y, r$X, binomial, r$beta, r$var, r$phi)$QIC
    index <- which.min(qic)
    select[1, index] <- select[1, index] + 1
  }
  colnames(select) <- c("  x1  ", "  x1 x2  ", "  x1 x3  ", "  x1 x2 x3  ", "  x1 x2 x3 x4  ")
  rownames(select) <- "time"
  return(select)
}

# variable(50, 100, "exchangeable")
