rm(list=ls())

U<-function(q){
  # s= if(sign(q)==1) sign(q) else 0 
  s = sapply(q, function(x) if(sign(x)==1) sign(x) else 0 )
  prob=c(0.3,0.2)
  U<- -sum(dnorm(q,0,1,log = TRUE)+dbinom(x=s,size = 1,prob = prob,log = TRUE))
  return(U)
}

# q<-c(1,2,-3,2,-4)
# u<-U(q)
# u
# U(q[1])+U(q[2])+U(q[3])+U(q[4])+U(q[5])

grad_U<-function(q){
  return(q)
}
# traceback()

res<- FIRSTDISCONTINUITY(c(0.0001235201,0.0001235),c(-1.6974357066,-1.697),0,U,0.025)

# FIRSTDISCONTINUITY(0.0001235,-1.697,0,U,0.025)
# t=epsilon-t0


FIRSTDISCONTINUITY <- function(q,p,t,U,epsilon){
  result<-matrix(NA,length(q),4)
  current_q<- q
  for (i in 1:length(q)) {
    q<-current_q[i]
    q_vector = c(q)
    num=0
    for (x in rep(epsilon/100,100)) {
      q_old = q_vector[length(q_vector)]
      q_new = q_old + p[i] * x
      q_vector=c(q_vector,q_new)
      num=num+1
      if ((num<100)& (q_old * q_new<=0)){
        break
      }
    }
    
    q_left_2 = q_vector[length(q_vector)-1]- p[i]*x # 为了计算方向
    q_left_1 = q_vector[length(q_vector)-1]# 为了计算方向
    q_right = q_vector[length(q_vector)]
    # print(q_vector)
    U_left_2 = U(q_left_2)
    U_left_1 = U(q_left_1)
    U_right = U(q_right)
    delta_U <- U_right - U_left_1
    tan_direction=(U_left_1 - U_left_2)/(epsilon/100)  # 方向的正切值
    if (num<100) {
      result[i,]= c( q_left_1, num*epsilon/100, delta_U,  tan_direction )
    }else{
      result[i,]=   rep(0,4)
    }
  }
  result<-data.frame(result)
  colnames(result)<-c('q','t_x','delta_U','tan')
  return (result)
}




# q0 ,current sample (discret parameter)
#  U ,potential function
# epsilon leapfrog step size
#  L  leapfrog steps 

new = RHMC(U = U,grad_U = grad_U,epsilon = 0.01,L=10,current_q = c(0.1,0.2))
new

RHMC<-function(U, grad_U, epsilon, L, current_q){
  q = current_q
  p = rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p
  for (i in 1:L) {
    # print('i')
    # print(i)
    p = p - epsilon * grad_U(q) / 2
    t0<-0
    
    discounity <- FIRSTDISCONTINUITY(q,p,t0,U,epsilon)
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
      
      for (j in 1:length(q)) {
        if ((K[j]/2) > delta_U[j]) {
          p[j] = sqrt(K[j]-2*delta_U[j]+p_parallel[j]^2)*sign(p[j])
        }else{
          # p_vetical = -p_vetical
          p[j] = -p[j]  # 反弹回去了
        }
      }

      
      # p = sqrt(p_vetical^2+p_parallel^2)
      
      q = q + (epsilon - t0) * p
    }else{
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

x0<-c(1,1)
bi<-10000
n<-30000
acc_hmc<-0

acc<-rep(0,n)
x<-matrix(NA,nrow = n,ncol = 2)
x[1,]= x0

for (i in 2:n) {
  old<-x[i-1,]
  # print('old')
  # print(old)
  new = RHMC(U = U,grad_U = grad_U,epsilon = 0.015,L=20,current_q = old)
  # print(new)
  if(new!=old && i>bi) {acc_hmc=acc_hmc+1}
  if(new!=old ){acc[i-1]=1} 
  x[i,]=new
  # print('\t')
  # print(i)
  if(i%%100==0){print(i)}
}
acc_hmc
length(x)
res=x[bi:n,]
s_res=sign(res)

hist(sign(res[,1]),100)
hist(sign(res[,2]),100)

library(cumstats)
plot(1:length(res[,1]),cummean(res[,1]),type='l')
plot(1:length(res[,2]),cummean(res[,2]),type='l')



new_res <- ifelse(res<0,0,1)
apply(new_res,2,mean)

plot(1:length(new_res[,1]),cummean(new_res[,1]),type='l')
plot(1:length(new_res[,2]),cummean(new_res[,2]),type='l')
tail(cummean(new_res[,1]))
tail(cummean(new_res[,2]))











