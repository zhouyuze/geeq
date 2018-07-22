generatedata <- function(beta,alpha,gamma,X,T,n)  {
  mean.vec <- exp(crossprod(t(X),beta))
  y <- matrix(0,nrow=n,ncol=T)
  y[,1] <- rnbinom(n,mu = mean.vec[1],size=mean.vec[1]/gamma)
  for (i in 1:n)  {
    for (t in 2:T)  {
      innovation.mean <- mean.vec[t] - alpha*(sqrt(mean.vec[t]*mean.vec[t-1]))
      I <- rnbinom(1,mu= innovation.mean,size= innovation.mean/gamma)
      first.shape <- alpha*sqrt(mean.vec[t]*mean.vec[t-1])/gamma
      second.shape <- mean.vec[t-1]/gamma - first.shape
      u <- rbeta(1,shape1 = first.shape,shape2=second.shape)
      a <- rbinom(1,size=y[i,t-1],prob=u)
      y[i,t] = a + I
    }
  }
  longform <- c(t(y))
  simdata <- data.frame(count = longform, time = rep(X[,2],times=n),subject=rep(c(1:n),each=T))
  return(simdata)
}

compare <- function(cluster.size) {
  X <- cbind(rep(1,5),c(-.5,-.25,0,.25,.5))
  group <- length(cluster.size)
  result <- data.frame(cluster.size, geeq=rep(0, group), qif=rep(0, group),
                       geepack=rep(0, group), geeM=rep(0, group))

  for (i in 1:group) {
    size = cluster.size[i]
    testdat <- generatedata(beta=c(1,.5),alpha=.2,gamma=.5,X=X,T=5,n=size)

    start <- Sys.time()
    gee(count~time, data=testdat, id=subject, family=poisson, corstr="ar1")
    end <- Sys.time()
    result[i, 2] <- as.double(end - start)

    start <- Sys.time()
    qif(count~time, data=testdat, id=subject, family=poisson, corstr="ar1")
    end <- Sys.time()
    result[i, 3] <- as.double(end - start)

    start <- Sys.time()
    geese(count~time, data=testdat, id=subject, family=poisson, corstr="ar1")
    end <- Sys.time()
    result[i, 4] <- as.double(end - start)

    start <- Sys.time()
    geem(count~time, data=testdat, id=subject, family=poisson, corstr="ar1")
    end <- Sys.time()
    result[i, 5] <- as.double(end - start)
  }
  result
}

cluster.size <- c(100, 200, 500, 1000, 2000, 5000)

compare(cluster.size)


