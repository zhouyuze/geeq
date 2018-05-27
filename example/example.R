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
  print(apply(y,2,mean))
  simdata <- data.frame(count = longform, time = rep(X[,2],times=n),subject=rep(c(1:n),each=T))
  return(simdata)
}

X <- cbind(rep(1,5),c(-.5,-.25,0,.25,.5))
testdat <- generatedata(beta=c(1,.5),alpha=.2,gamma=.5,X=X,T=5,n=200)

y <- testdat$count
X <- cbind(rep(1, 1000), testdat$time)
clusterSize <- rep(5, 200)
family <- poisson()
qif_c(y, X, rep(0, length(y)), clusterSize, family, "ar1", c(1.0205254,0.4931009))
