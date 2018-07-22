library(bindata)
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

correlation_select <- function(n, repeat_time) {
  qic <- c(0, 0, 0)
  select <- matrix(0, ncol=3, nrow=1)
  for (i in 1:repeat_time) {
    data <- generatedata(n, 3)
    qic[1] <- gee(y~x1+x2, data, id, family=binomial, corstr="independence")$QIC
    qic[2] <- gee(y~x1+x2, data, id, family=binomial, corstr="exchangeable")$QIC
    qic[3] <- gee(y~x1+x2, data, id, family=binomial, corstr="ar1")$QIC
    index <- which.min(qic)
    select[1, index] <- select[1, index] + 1
  }
  colnames(select) <- c("independence", "exchangeable", "ar1")
  rownames(select) <- "time"
  return(select)
}

variable_select <- function(n, repeat_time, cor_str) {
  qic <- rep(0, 5)
  select <- matrix(0, ncol=5, nrow=1)
  for (i in 1:repeat_time) {
    data <- generatedata(n, 3)
    qic[1] <- gee(y~x1, data, id, family=binomial, corstr=cor_str)$QIC
    qic[2] <- gee(y~x1+x2, data, id, family=binomial, corstr=cor_str)$QIC
    qic[3] <- gee(y~x1+x3, data, id, family=binomial, corstr=cor_str)$QIC
    qic[4] <- gee(y~x1+x2+x3, data, id, family=binomial, corstr=cor_str)$QIC
    qic[5] <- gee(y~x1+x2+x3+x4, data, id, family=binomial, corstr=cor_str)$QIC
    index <- which.min(qic)
    select[1, index] <- select[1, index] + 1
  }
  colnames(select) <- c("x1", "x1 x2", "x1 x3", "x1 x2 x3", "x1 x2 x3 x4")
  rownames(select) <- "time"
  return(select)
}

correlation_select(50, 1000)
variable_select(50, 1000, "independence")
variable_select(50, 1000, "exchangeable")
variable_select(50, 1000, "ar1")

correlation_select(100, 1000)
variable_select(100, 1000, "independence")
variable_select(100, 1000, "exchangeable")
variable_select(100, 1000, "ar1")
