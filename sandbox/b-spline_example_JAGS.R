# Example of b-splines in JAGS from:
#  https://www4.stat.ncsu.edu/~reich/st590/code/Moto.html



# Libraries ---------------------------------------------------------------

library(splines)
library(rjags)
library(MASS)


# Data --------------------------------------------------------------------

Y <- mcycle$accel
X <- mcycle$times

Y <- (Y-mean(Y))/sd(Y)
X <- X/max(X)

n <- length(Y)


# Basis functions ---------------------------------------------------------

J <- 6
B <- bs(X,J)
dim(B)
matplot(X,B,type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)


# Example posterior draws -------------------------------------------------

g <- matrix(0,n,10)

for(j in 1:10){
  beta  <- rnorm(J, 0, 0.1)
  g[,j] <- B%*%beta
}

matplot(X,g,lwd=2,type="l",cex.lab=1.5,cex.axis=1.5)


# JAGS model --------------------------------------------------------------

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

# dat    <- list(Y=Y,n=n,B=B,J=J)
# init   <- list(mu=mean(Y),beta=rep(0,J),taue=var(Y))
# model <- jags.model(textConnection(Moto_model),inits=init,data = dat, n.chains=2)
# update(model, 10000, progress.bar="none")
# samp   <- coda.samples(model, 
#                        variable.names=c("mean", "mu"), 
#                        n.iter=20000, progress.bar="none")
# sum <- summary(samp)
# q <- sum$quantiles
# 
# plot(X,Y,xlab="time",ylab="Acceleration",cex.lab=1.5,cex.axis=1.5)
# 
# lines(X,q[,1],col=2,lty=2)
# lines(X,q[,3],col=2,lty=1)
# lines(X,q[,5],col=2,lty=2)
# legend("bottomright",c("Median","95% interval"),lty=1:2,col=2,bg=gray(1),inset=0.05,cex=1.5)


# Example with 26 week sin wave -------------------------------------------
dosin <- function(nu, phi, t, period = 26){
  return( (0 + nu * sin( (2*pi*t)/period + phi ) ) )
}

mywave <- dosin(0.2, 13, 1:26, 26) + rnorm(26, 0, 0.1)
plot(mywave)

times <- 1:26
knots <- 6
bases <- bs(times, knots)
matplot(times,bases,type="l",xlab="Time",ylab="Basis function, Bj(X)")

sin_model <- "model{

  # Likelihood
  for(i in 1:n){
    y[i] ~ dnorm(mean[i], taue)
    mean[i] <- mu + inprod(B[i,], beta[])
  }

  # Random effects, b-splines
  for(i in 1:J){
    beta[i] ~ dnorm(0, taub)
  }

  # Priors
  mu   ~ dnorm(0,0.01)
  taue ~ dgamma(0.1,0.1)
  taub ~ dgamma(0.1,0.1)

}"

dat <- list(y = mywave, n=length(mywave), B=bases, J=knots)
init   <- list(mu=mean(mywave),beta=rep(0,knots),taue=var(mywave))
model <- jags.model(textConnection(sin_model),inits=init,data = dat, n.chains=2)
update(model, 10000, progress.bar="none")
samp   <- coda.samples(model,
                       variable.names=c("mean"),
                       n.iter=20000, progress.bar="none")

summ <- summary(samp)
q <- summ$quantiles

plot(times,mywave,xlab="time",ylab="beta")
lines(times,q[,1],col=2,lty=2)
lines(times,q[,3],col=2,lty=1)
lines(times,q[,5],col=2,lty=2)
legend("bottomright",c("Median","95% interval"),lty=1:2,col=2,bg=gray(1),inset=0.05)

