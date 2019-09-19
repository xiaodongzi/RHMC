rm(list=ls())

U<-function(q){
  s= if(sign(q)==1) sign(q) else 0 
  U<- -dnorm(q,0,3,log = TRUE)- dbinom(x=s,size = 1,prob = 0.3,log = TRUE)
  return(U)
}


grad_U<-function(q){
  return(q)
}


# traceback()

# FIRSTDISCONTINUITY(0.0001235201,-1.6974357066,0,U,0.025)
# FIRSTDISCONTINUITY(1,-2,0,U,0.2)

# t=epsilon-t0

FIRSTDISCONTINUITY<- function(q,p,t,U,epsilon){
  q_vector = c(q)
  num=0
  iter = 50
  for (x in rep(epsilon/iter,iter)) {
    q_old = q_vector[length(q_vector)]
    q_new = q_old + p * x
    q_vector=c(q_vector,q_new)
    num=num+1
    if ((num<iter)& (q_old * q_new<=0)){
      # print(c(q_vector,p))
      break
    }
  }
  
  q_left_2 = q_vector[length(q_vector)-1]- p*x # 为了计算方向
  q_left_1 = q_vector[length(q_vector)-1]# 为了计算方向
  q_right = q_vector[length(q_vector)]
  # print(q_vector)
  U_left_2 = U(q_left_2)
  U_left_1 = U(q_left_1)
  U_right = U(q_right)
  delta_U <- U_right - U_left_1
  tan_direction=(U_left_1 - U_left_2)/(epsilon/iter)  # 方向的正切值
  if (num<iter) {
    return(list(q=q_left_1, t_x = num*epsilon/iter, delta_U= delta_U, tan = tan_direction ))
  }else{
    return(0)    
  }
  
}



# q0 ,current sample (discret parameter)
#  U ,potential function
# epsilon leapfrog step size
#  L  leapfrog steps 
# new = RHMC(U = U,grad_U = grad_U,epsilon = 0.25,L=25,current_q = old)
# new

RHMC<-function(U, grad_U, epsilon, L, current_q){
  q = current_q
  p = rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p
  for (i in 1:L) {
    # print('i')
    # print(i)
    p = p - epsilon * grad_U(q) / 2
    t0<-0
    discounity<-0
    if (abs(q)<epsilon){
    discounity <- FIRSTDISCONTINUITY(q,p,t0,U,epsilon)}
    # print(discounity)
    if (discounity!=0) {
      q = discounity$q
      # print('q')
      # print(q)
      t0 = discounity$t_x
      delta_U= discounity$delta_U
      
      tan = discounity$tan  #  如果能量往下走，正切可能是负值
      p_parallel <- p* abs(tan/ sqrt(1+tan^2))  # 平行=p* sin
      # print('p_parallel')
      # print(p_parallel)
      p_vetical <- p /sqrt(1+tan^2)
      # print('p_vetical')
      # print(p_vetical)
      K<- p_vetical^2
      # print('K')
      # print(K)
      if ((K/2) > delta_U) {
        p_vetical<- sqrt(K-2*delta_U)
        p = sqrt(p_vetical^2+p_parallel^2)*sign(p)
      }else{
        # p_vetical = -p_vetical
        p = -p  # 反弹回去了
      }
      
      # p = sqrt(p_vetical^2+p_parallel^2)
      
      q = q + (epsilon - t0) * p
    }
    else{
      q = q + epsilon  * p
    }
    
    # print('qqqq')
    # print(q)
    p = p - epsilon * grad_U(q) / 2
  }
  p = -p
  
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (q)  # accept
  }
  else
  {
    return (current_q)  # reject
  }
}

#  接下来开始仿真模拟 

x0<- 1
bi<-1000
n<-20000
acc_hmc<-0

acc<-rep(0,n)
x<-rep(NA,n)
x[1]=x0
for (i in 2:n) {
  old<-ifelse((i==1),x0,x[i-1])
  # print('old')
  # print(old)
  
  new=RHMC(U = U,grad_U = grad_U,epsilon = runif(1,0.1,0.9),L=runif(1,10,20),current_q = old)
  if(new!=old && i>bi) {acc_hmc=acc_hmc+1}
  if(new!=old){acc[i-1]=1} 
  x[i]=new
  # print('\t')
  # print(i)
  if(i%%100==0){print(i)}
}
acc_hmc
length(x)
res=x[bi:n]
hist(sign(res),100)
library(cumstats)
plot(1:length(res),cummean(res),type='l')


new_res=c()
for (x in res){
  if (x<=0) {
    new_res=c(new_res,0)
  }else{
    new_res=c(new_res,1)
  }
}
plot(1:length(new_res),cummean(new_res),type='l')

tail(cummean(new_res))

mean(new_res)




