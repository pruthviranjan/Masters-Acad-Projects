library(mlmRev)
library(ggplot2)

set.seed(7052)
#load dataset
early <- data.frame(Early)
#google high cog score > 110
early$high <- ifelse(early$cog >= 110, 1, 0)

#cummulative distribution of cog scores of infants age 1.5 yrs
inf_1.5 <- early[early$age == 1.5,]
ggplot(data=inf_1.5, aes(x = cog))+geom_histogram(bins = 8)

cog15 <- ecdf(inf_1.5$cog)

Alpha <- 0.05
n <- length(inf_1.5$cog)
eps <- sqrt((1/(2*n))*log(2/Alpha))
grid <- seq(60,150,length.out = 1000)
plot(cog15)
lines(grid,pmin(cog15(grid) + eps,1))
lines(grid,pmax(cog15(grid) - eps,0))

treatment <- inf_1.5$cog[inf_1.5$trt == "Y"]
control<- inf_1.5$cog[inf_1.5$trt == "N"]

#mle and asymptotic dist

mean_treat.hat<-mean(treatment)
mean_control.hat<-mean(control)
sd_treat.hat<-sd(treatment)
sd_control.hat<-sd(control)

mle <- mean_treat.hat - mean_control.hat
print(mle)
## [1] 14.40421
asym.sd <- sqrt((var(treatment) + var(control))/58)
ci <- c(mle - 2*asym.sd, mle + 2*asym.sd)
print(ci)
## [1]  9.192859 19.615570
#bootstrap and CI

B <- 1000
theta.par <- rep(0,B)
theta.nonpar <- rep(0,B)

for(i in 1:B){
  x1 <- rnorm(58, mean_treat.hat, sd_treat.hat)
  x2 <- rnorm(58, mean_control.hat, sd_control.hat)
  theta.par[i] <- mean(x1) - mean(x2)
  
  x1b <- sample(treatment,size = 58, replace = T)
  x2b <- sample(control, size = 58, replace = T)
  theta.nonpar[i] <- mean(x1b) - mean(x2b)
}

sd.par <- sd(theta.par)
print(sd.par)
## [1] 2.615894
sd.nonpar <- sd(theta.nonpar)
print(sd.nonpar)
## [1] 2.698259
normal.ci<-c(mle-2*sd.nonpar, mle+2*sd.nonpar)
print(normal.ci)
## [1]  9.007697 19.800732
pivotol.ci<-c(2*mle-quantile(theta.nonpar,0.975), 2*mle-quantile(theta.nonpar,0.025))
print(pivotol.ci)
##     97.5%      2.5% 
##  9.169636 19.842912
quantile.ci<-quantile(theta.nonpar, c(0.025, 0.975))
print(quantile.ci)
##      2.5%     97.5% 
##  8.965517 19.638793
par.ci <- c(mle - 2*sd.par, mle + 2*sd.par)
print(par.ci)
## [1]  9.172427 19.636002
#hypothesis test for treatment vs control groups of infants aged 1.5
w.stat <- mle/asym.sd
pvalue <- 2*pnorm(-abs(w.stat))
print(pvalue)
## [1] 3.238825e-08
#at 95% confidence, we reject null
ggplot(inf_1.5, aes(cog, fill = trt)) + geom_density(alpha = 0.2)

#Bayesian
quantile.5 <- list(p=0.5, x=0.25)
quantile.99 <- list(p=0.99999,x=0.40)
quantile.0 <- list(p=0.00001,x=0.05)


library(LearnBayes)

prior.1<-beta.select(quantile.5,quantile.99)
prior.2<-beta.select(quantile.5, quantile.0)

max.error<-10000000
for(i in seq(from=min(prior.1[1],prior.2[1]),to=max(prior.1[1],prior.2[1]),length.out = 100))
{
  for(j in seq(from=min(prior.1[2],prior.2[2]),to=max(prior.1[2],prior.2[2]),length.out = 100))
  {
    prior.calc.1=qbeta(p = 0.5,shape1 = i,shape2 = j)
    prior.calc.2=qbeta(p = 0.999,shape1 = i,shape2 = j)
    prior.calc.3=qbeta(p = 0.001,shape1 = i,shape2 = j)
    prior.error<-abs(prior.calc.1-0.15)+abs(prior.calc.2-0.25)+abs(prior.calc.3-0.05)
    if(prior.error<max.error){
      max.error<-prior.error
      alpha <-i
      beta <-j
    }
  }
}
plot(ylab = "Density", xlab ="theta" ,seq(from = 0,to = 1,length.out = 1000),dbeta(x = seq(from = 0,to = 1,length.out = 1000),shape1 = alpha, shape2 = beta))

#posterior
library(dplyr)
subset<-sample_n(tbl = inf_1.5,size = 100)
s <-sum(subset$cog>=110)
n <-nrow(subset)

posterior.alpha <- alpha + s
posterior.beta <- beta + n - s


posterior <- dbeta(x = seq(0.005, 1, length = 5000), posterior.alpha, posterior.beta)
prior<- dbeta(x = seq(0.005, 1, length = 5000),shape1 = alpha, shape2 = beta)
plot(seq(0.005, 1, length = 5000), posterior,col="red")
lines(seq(0.005, 1, length = 5000), prior, col = "blue")
legend(x = 0.8,y = 12,legend = c("Prior","Posterior"),lwd=c(1,1),col=c("red","blue"))

#posterior mean
posterior.mean <- posterior.alpha/(posterior.alpha + posterior.beta)
print(posterior.mean)
## [1] 0.2279329
posterior.var <- (posterior.alpha*posterior.beta)/((posterior.alpha + posterior.beta)^2)*(posterior.alpha + posterior.beta + 1)
print(posterior.var)
## [1] 36.23923
