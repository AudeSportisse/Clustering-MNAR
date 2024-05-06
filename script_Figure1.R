library("MASS")
library(mvtnorm)

K=3
mu=c(-5,0,5)
var=c(1,2,3)
pi=c(0.3,0.3,0.4)
alpha=c(2,0,1)

rho <- function(alpha){
  return(1/(1+exp(-alpha)))
}


proba_c_given_y <- function(x){
  proba_z_given_y <- numeric(K)
  for (k in 1:K){
    proba_z_given_y[k] <- dnorm(x,mean=mu[k],sd=sqrt(var[k]))*pi[k]
  }
  proba_z_given_y <- proba_z_given_y/sum(proba_z_given_y)
  
  proba_c_given_y <- numeric(K)
  for (k in 1:K){
    proba_c_given_z <- rho(alpha[k])
    proba_c_given_y[k] <- proba_c_given_z*proba_z_given_y[k]
  }
  return(sum(proba_c_given_y))
}

y <- seq(from=-10, to=10, length.out=100)
plot(y,sapply(y,proba_c_given_y),type = "l", lty = 1)


##GGplot

library(ggplot2)
library(latex2exp)
df_plot_main = data.frame(abs=y,ord=sapply(y,proba_c_given_y))
ggplot(df_plot_main, aes(x=abs, y=ord)) + geom_line() +
  ylab(TeX("$P(c=1|y)$")) + xlab("y")
                    