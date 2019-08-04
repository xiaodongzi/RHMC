#  模拟抽样t 分布
rm(list = ls())
x0<-0
bi<-10000
n<-110000
acc_rand<-0
x<-rep(NA,n)
nu<-2
#  restore x
for (i in 1:n) {
  old<-ifelse((i==1),x0,x[i-1])
  new<-rnorm(1,old,nu)
  alpha<-min(1,dt(new,5)/dt(old,5))
  accept<-rbinom(1,1,alpha)
  if(accept){x[i]=new}
  if(accept && i>bi)(acc_rand=acc_rand+1)
  if(!accept){x[i]=old}
}
result<-x[10001:110000]
acc_rand



summary(result)

real=rt(n-bi,5)
summary(real)



#  hmc 
U<-function(q){
  q_energy=3*log(1+ q^2/5)+0.9686196
  return(q_energy)
}


grad_u<-function(q){
  # grad=-6*q*exp(log(1+q^2/5)*(-4))
  grad= 6*q/(5+q^2)
  return(grad)
}

HMC = function (U, grad_U, epsilon, L, current_q)
{
  q = current_q
  p = rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p
  
  # Make a half step for momentum at the beginning
  
  p = p - epsilon * grad_U(q) / 2
  
  # Alternate full steps for position and momentum
  
  for (i in 1:L)
  {
    # Make a full step for the position
    
    q = q + epsilon * p
    
    # Make a full step for the momentum, except at end of trajectory
    
    if (i!=L) p = p - epsilon * grad_U(q)
  }
  
  # Make a half step for momentum at the end.
  
  p = p - epsilon * grad_U(q) / 2
  
  # Negate momentum at end of trajectory to make the proposal symmetric
  
  p = -p
  
  # Evaluate potential and kinetic energies at start and end of trajectory
  
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
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

x0<-0
bi<-10000
n<-110000
acc_hmc<-0

acc<-rep(0,n)
x<-rep(NA,n)
x[1]=x0
for (i in 2:n) {
  old<-ifelse((i==1),x0,x[i-1])
  new=HMC(U = U,grad_U = grad_u,epsilon = 0.02,L=25,current_q = old)
  if(new!=old && i>bi){acc_hmc=acc_hmc+1}
  if(new!=old){acc[i-1]=1} 
  x[i]=new
  if(i%%1000==0){print(i)}
}
hmc<-x[10001:110000]
hist(hmc,100)
summary(hmc)
summary(real)
summary(result)
acc_hmc
hist(real,100)

library(cumstats)
plot(1:length(cummean(result)),cummean(result),type='l')


plot(1:100000,hmc,'l')


test=seq(-5,5,0.01)
test_1=dnorm(test)
plot(test,test_1,'l')
test_2<--log(test_1)
plot(test,test_2,'l')



