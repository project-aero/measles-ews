t <- 100
mu <- log(0.00003)
eps <- rnorm(t, 0, 1)
f <- 0.05
# plot(exp(mu*(1+f*eps)), type = "l")


mean(rbinom(100000, size = 1, prob = 0.5))
median(rbinom(100000, size = 1, prob = 0.5))



ar1log <- function(x0, gens, ar, p){
  x <- numeric(gens)
  x[1] <- x0
  for(i in 2:gens){
    x[i] <- x[i-1]*ar + p*rnorm(1,0,1)
  }
  return(x)
}

xout <- matrix(nrow = 100, ncol = 1000)
for(i in 1:100){
  xout[,i] <- ar1log(log(0.00003), gens = 100, ar = 1, p = 0.05)
}

matplot(exp(xout), type = "l", col = "grey", lty = 1)
