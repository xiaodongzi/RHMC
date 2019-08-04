#  HMC  regression 
rm(list = ls())
#  data generate process 
n=2000
p<-2
epsilon<-rnorm(n,0,2)
x<-matrix(runif(n*p),n)
x<-cbind(1,x)
beta_dgp<-matrix(c(-3,1,3))
y<- x %*% beta_dgp + epsilon

a<-lm(y~x)
summary(a)
library(mvtnorm)



U_Posterior<-function(q,y,x){
  beta<-matrix(q[1:3])
  logsigma2<-q[4]
  logprior_beta<- dmvnorm(t(beta),#  是行向量
                         mean = matrix(0,1,3),
                         sigma = diag(3),
                         log=T)
  
  logprior_logsigma2<- dnorm(logsigma2,0,1,log = T)
  
  loglikelihood<-sum(dnorm(x = y,
                           mean = x %*% beta,
                           sd = sqrt(exp(logsigma2)),
                           log = T))
  logposterior<-logprior_beta+logprior_logsigma2+loglikelihood
  
  U<- -logposterior
  # U<--1*(log(sigma2)*(-n/2+4)- 1/2*t(beta)%*%beta-sigma2/2-
  # 1/(2*sigma2)*sum(eplison^2))
  return(U)
}


grad_u<-function(q,y,x){
  beta<-matrix(q[1:3])
  sigma2<-exp(q[4])
  eplison<- y-x%*%beta
  n=length(y)
  logU_beta0<- -1*(-beta[1]+ 1/sigma2*sum(eplison*x[,1]))
  logU_beta1<- -1*(-beta[2]+ 1/sigma2*sum(eplison*x[,2]))
  logU_beta2<- -1*(-beta[3]+ 1/sigma2*sum(eplison*x[,3]))
  # logU_sigma2<- -1*((4-n/2)/sigma2-1/2+1/(2*sigma2^2)*sum(eplison^2))
  logU_sigma2<- -1*(-log(sigma2)-n/2+1/(2*sigma2)*sum(eplison^2))
  
  grad<-matrix(c(logU_beta0,logU_beta1,logU_beta2,logU_sigma2))
  return(grad)
}



HMC = function (U, grad_U, epsilon, L, current_q,y,x){
  q = current_q
  p = rnorm(length(q),0,2)  # independent standard normal variates
  current_p = p
  
  # Make a half step for momentum at the beginning
  
  p = p - epsilon * grad_U(q,y,x) / 2
  
  # Alternate full steps for position and momentum
  
  for (i in 1:L)
  {
    # Make a full step for the position
    
    q = q + epsilon * p
    
    # Make a full step for the momentum, except at end of trajectory
    
    if (i!=L) p = p - epsilon * grad_U(q,y,x)
  }
  
  # Make a half step for momentum at the end.
  
  p = p - epsilon * grad_U(q,y,x) / 2
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q,y,x)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q,y,x)
  proposed_K = sum(p^2) / 2
  # print(current_U)
  # print(proposed_U)
  # print(current_K)
  # print(proposed_K)
  # print(exp(current_U-proposed_U+current_K-proposed_K))
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (q)  # accept
  }
  else
  {
    return (current_q)  # reject
  }
}


beta_init<-c(1,1,1)
sigma_init<-1
niter<-21000
bi<-10000
para<-matrix(NA,niter,4)
para[1,]<-c(beta_init,sigma_init)
acc_hmc<-0

for (i in 2:niter) {
  if(i%%100==0){print(i)}
  old<- para[i-1,]
  new=HMC(U = U_Posterior,grad_U = grad_u,epsilon = 0.002,L = 50,current_q = old,y,x)
  if(new!=old && i>bi){acc_hmc=acc_hmc+1}
  para[i,]=new
}

result=para[bi:niter,]
apply(result,2,mean)
acc_hmc
hist(result[,1],100)
hist(result[,2],100)
hist(result[,3],100)
hist(result[,4],100)

library(cumstats)
plot(1:length(cummean(result[,1])),cummean(result[,1]),type='l')  
plot(1:length(cummean(result[,2])),cummean(result[,2]),type='l')   
plot(1:length(cummean(result[,3])),cummean(result[,3]),type='l')   
plot(1:length(cummean(result[,1])),cummean(result[,4]),type='l')   







