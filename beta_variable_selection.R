#  spike and slab regression


rm(list=ls())
library(truncnorm)
library(mvtnorm)
n=1000
noise=rnorm(n,0,1)
col=5
x=rnorm(n*col,0,1)
X=matrix(x,n)
beta<-c(-2,2,3,4,-3)
Y= X %*% beta + noise
extra_x=rnorm(n,0,1)
X=cbind(X,extra_x)

#   use ols finding extra_x is not significant
exp=lm(Y~X)
summary(exp)
# 
col=col+1
beta_0=c(1,0,1,1,0,2)
spike_p_0=rep(0.5,col)
s_0=c(1,0,1,1,0,1)


LogPosterior<-function(spike_p,s,beta,X,Y){
  temp=which(s!=0)
  # print(cat('temp',temp))
  X=X[,temp]
  beta=beta[temp]
  col=length(beta)
  log_spike_p_prior<-0 # all is runif
  # log_spike_p_prior<-sum(log(dtruncnorm(spike_p*10,a=0,b=10,mean=5,sd=4)))
  log_spike_s_prior<- sum(dbinom(s,size=1,prob = spike_p,log = TRUE))
  
  Log_Likelihood<- sum(dnorm(x=Y,
                             mean = X %*% matrix(beta),
                             sd = 1,
                             log=TRUE))
  
  Log_beta_prior<-dmvnorm(x= beta,
                          mean = matrix(0,1,col),
                          sigma = diag(col),
                          log = TRUE)
  
  Log_posterior<-Log_beta_prior + Log_Likelihood + log_spike_s_prior + log_spike_p_prior
  return(Log_posterior)
}

nIter=50000
spike_p<-matrix(NA,nIter,col)
acc<-matrix(NA,nIter,1)
spike_p[1,]= spike_p_0
beta<-matrix(NA,nIter,col)
beta[1,]=beta_0
s<-matrix(NA,nIter,col)
s[1,]=s_0

for (i in 2:nIter){
  
  if(i%%100==0){print(i)}
  p.current = spike_p[i-1,]
  beta.current = beta[i-1,]
  s.current = s[i-1,]
  
  # the metropolis for spike_p
  # p.prop = runif(col)
  p.prop = rbeta(n= col ,1,1)
  
  LogPosterior_p.prop <- LogPosterior(spike_p = p.prop,
                                      s =s.current,
                                      beta = beta.current,
                                      X = X,Y= Y )
  Log_hasting_p.prop = sum(log(dtruncnorm(p.prop ,a = 0,b=1,mean = 0.5,sd = 0.2)))
  
  LogPosterior_p.current<-LogPosterior(spike_p = p.current,
                                       s=s.current,
                                       beta=beta.current,
                                       X = X, Y = Y)
  Log_hasting_p.current = sum(log(dtruncnorm(p.current,a = 0,b=1,mean = 0.5,sd = 0.2)))
  
  # logratio<-LogPosterior_p.prop - LogPosterior_p.current
  logratio<-LogPosterior_p.prop - LogPosterior_p.current-Log_hasting_p.prop + Log_hasting_p.current
  
  acc.prop.p<- min(exp(logratio),1)
  # print(cat('logratio',logratio))
  if(runif(1)< acc.prop.p){
    p.current= p.prop
  }
  spike_p[i,]=p.current
  
  #  the metropolis for s 
  s.prop<- rbinom(n = 6,size = 1,prob = p.current)
  LogPosterior_s.prop=LogPosterior(spike_p = p.current,
                                   s=s.prop,
                                   beta = beta.current,
                                   X = X,Y = Y )
  Log_hasting_s.prop = dbinom(x = s.prop,size = 1,prob = p.current,log = TRUE)
  
  LogPosterior_s.current<-LogPosterior(spike_p = p.current,
                                       s=s.current,
                                       beta = beta.current,
                                       X = X,Y = Y )
  Log_hasting_s.current<-dbinom(x=s.current,size = 1,prob=p.current,log = TRUE)
  
  logratio=LogPosterior_s.prop - LogPosterior_s.current + Log_hasting_s.current-Log_hasting_s.prop
  
  acc.prop.s<-min(exp(logratio),1)
  if(runif(1)< acc.prop.s){
    s.current= s.prop
  }
  s[i,]=s.current
  
  #  the metropolis for beta
  beta.prop<-rmvnorm(n=1,
                     mean = matrix(beta.current,1),
                     sigma = diag(col)
  )
  
  
  # beta.prop.zero_index=which(s.current==0)
  # beta.prop[beta.prop.zero_index]=NA
  # beta.prop=as.vector(beta.prop)
  
  
  LogPosterior_beta.prop<-LogPosterior(spike_p = p.current,
                                       s=s.current,
                                       beta = beta.prop,
                                       X=X,Y=Y)
  
  LogPosterior_beta.current<-LogPosterior(spike_p = p.current,
                                          s=s.current,
                                          beta= beta.current,
                                          X=X,Y=Y)
  logratio<-LogPosterior_beta.prop - LogPosterior_beta.current
  
  acc.prop.beta<- min(exp(logratio),1)
  
  if(runif(1)<acc.prop.beta){
    beta.current<-beta.prop
  }
  beta[i,]= beta.current
  
  
}

# 主要问题是现在要处理，当s=0 时候beta 不是全部都要用上的。
beta_mean=apply(beta,2,mean)
print('beta_mean')
beta_mean
s_mean=apply(s,2,mean)
print('s_mean')
s_mean
p_mean=apply(spike_p,2,mean)
print('p_mean')
p_mean



hist(beta[,1],100,xlim = c(-5,5))
hist(beta[,2],100,xlim = c(-5,5))
hist(beta[,3],100,xlim = c(-5,5))
hist(beta[,4],100,xlim = c(-5,5))
hist(beta[,5],100,xlim = c(-5,5))
hist(beta[,6],100,xlim = c(-5,5))
#
hist(spike_p[,1],100,xlim = c(0,1))
hist(spike_p[,2],100,xlim = c(0,1))
hist(spike_p[,3],100,xlim = c(0,1))
hist(spike_p[,4],100,xlim = c(0,1))
hist(spike_p[,5],100,xlim = c(0,1))
hist(spike_p[,6],100,xlim = c(0,1))
#
hist(s[,1],100,xlim = c(0,1))
hist(s[,2],100,xlim = c(0,1))
hist(s[,3],100,xlim = c(0,1))
hist(s[,4],100,xlim = c(0,1))
hist(s[,5],100,xlim = c(0,1))
hist(s[,6],100,xlim = c(0,1))

# after_beta=beta[50000:100000,]
# after_s=s[50000:100000,]

library(cumstats)

plot(1:length(cummean(beta[10000:50000])),cummean(beta[10000:50000]),type = 'l')

library(coda)


temp<-apply(s,1,paste,collapse = ' ' )
table(temp)


beta_res=beta[temp =="1 1 1 1 1 0",]


chain=mcmc(beta_res)
plot(chain)

library(cumstats)
plot(1:length(cummean(beta_res[,1])),cummean(beta_res[,1]),type='l')   

plot(1:length(cummean(beta_res[,1])),cummean(beta_res[,2]),type='l')   
plot(1:length(cummean(beta_res[,1])),cummean(beta_res[,3]),type='l')   

plot(1:length(cummean(beta_res[,1])),cummean(beta_res[,4]),type='l')   
plot(1:length(cummean(beta_res[,1])),cummean(beta_res[,5]),type='l')   
plot(1:length(cummean(beta_res[,1])),cummean(beta_res[,6]),type='l')   





