
# Example of b-splines in JAGS from:
#  https://www4.stat.ncsu.edu/~reich/st590/code/Moto.html

library(splines)
library(rjags)
library(MASS)

Y <- mcycle$accel
X <- mcycle$times

Y <- (Y-mean(Y))/sd(Y)
X <- X/max(X)

n <- length(Y)
n

J <- 6
B <- bs(X,J)
dim(B)
matplot(X,B,type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)


g <- matrix(0,n,10)

for(j in 1:10){
  beta  <- rnorm(J, 0, 0.1)
  g[,j] <- B%*%beta
}

matplot(X,g,lwd=2,type="l",cex.lab=1.5,cex.axis=1.5)

Moto_model <- "model{

   # Likelihood
   for(i in 1:n){
      Y[i]    ~ dnorm(mean[i],taue)
      mean[i] <- mu + inprod(B[i,],beta[])
   }

   # Random effects
   for(i in 1:J){
    beta[i] ~ dnorm(0,taub)
   }

   # Priors
   mu   ~ dnorm(0,0.01)
   taue ~ dgamma(0.1,0.1)
   taub ~ dgamma(0.1,0.1)

  }"

dat    <- list(Y=Y,n=n,B=B,J=J)
init   <- list(mu=mean(Y),beta=rep(0,J),taue=var(Y))
model <- jags.model(textConnection(Moto_model),inits=init,data = dat, n.chains=2)
update(model, 10000, progress.bar="none")
samp   <- coda.samples(model, 
                       variable.names=c("mean", "mu"), 
                       n.iter=20000, progress.bar="none")
sum <- summary(samp)
q <- sum$quantiles

plot(X,Y,xlab="time",ylab="Acceleration",cex.lab=1.5,cex.axis=1.5)

lines(X,q[,1],col=2,lty=2)
lines(X,q[,3],col=2,lty=1)
lines(X,q[,5],col=2,lty=2)

legend("bottomright",c("Median","95% interval"),lty=1:2,col=2,bg=gray(1),inset=0.05,cex=1.5)

