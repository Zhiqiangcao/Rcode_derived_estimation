###Bias from multivariate, univariate and derived three methods of estimation based on Cox regression
###this program can reproduce simulations results in  Table S8 of supplementary material
#date: 2024.9.1
#author: Zhiqiang Cao

##generate data
rm(list=ls())
library(survival) 
library(mvtnorm)
library(MASS)

###############data generation################################
gend = function(n,c,b,mut,covt,cove,beta_q){
  b1 = b[1]; b2 = b[2]; b3 = b[3]
  #true variables
  datat = rmvnorm(n=n, mean=mut, sigma=covt)
  t1 = datat[,1]; t2 = datat[,2]; t3 = datat[,3] 
  u = runif(n)
  t = -log(u)/exp(b1*t1+b2*t2+b3*t3) #baseline hazard function lambda_0(t)=1
  cen = runif(n, min=0, max=c)
  delta = 1*(t<=cen)    
  time = delta*t+(1-delta)*cen  #observed time
  #error variables
  datae = rmvnorm(n=n, mean=c(0,0,0), sigma=cove) 
  e1 = datae[,1]; e2 = datae[,2]; e3 = datae[,3]
  #additive error model
  x1 = beta_q[1]*t1+e1; x2 = beta_q[2]*t2+e2
  x3 = beta_q[3]*t3+e3; 
  datat = data.frame(
    t = time,
    status = delta,
    t1 = t1,
    t2 = t2,
    t3 = t3,
    x1 = x1,
    x2 = x2,
    x3 = x3
  )
  return(datat)
}

###parameter settings
N = 1000  #simulation times
n = 400   #sample size

scenario_set = 1:2  
for(scenario in scenario_set){
  if(scenario == 1){  #for male
    trube = c(0.172,0.228,0.187) 
    mut = c(9.356683,8.654661,11.937024) 
    c = 0.005  #to produce 50% censoring rate
    beta_q = c(0.3667,0.4148,0.2762)
    sigma_t1 = sqrt(1.05927134)
    sigma_t2 = sqrt(0.676466609)
    sigma_t3 = sqrt(1.13869850)
    sigma_e1 = sqrt(0.8575957)
    sigma_e2 = sqrt(0.8836357)
    sigma_e3 = sqrt(0.9131195)
  }else if(scenario == 2){ #for female
    trube = c(0.160,0.112,0.074)
    mux = c(9.245602,7.510393,11.167943)
    c = 0.068  #to produce about 50% censoring rate
    beta_q = c(0.2980,0.2535,0.2515)
    sigma_t1 = sqrt(1.07312719)
    sigma_t2 = sqrt(1.07539124)
    sigma_t3 = sqrt(1.17795383)
    sigma_e1 = sqrt(0.9046743)
    sigma_e2 = sqrt(0.9308771)
    sigma_e3 = sqrt(0.9254783)
  }
  
  rho_t_set = seq(0.1,0.4,by=0.1)
  len_rhot = length(rho_t_set)
  rho_e_set = c(0.7,0.8,0.9,0.95)
  len_rhoe = length(rho_e_set)
  b1e_n = b2e_n = b3e_n = b1e_un = b2e_un = b3e_un = b1e_d = b2e_d = b3e_d = matrix(0,len_rhoe,len_rhot)
  
  for(k in 1:len_rhot){
    rho_t = rho_t_set[k]
    for(j in 1:len_rhoe){
      rho_e = rho_e_set[j]
      covt = matrix(c(sigma_t1^2,rho_t*sigma_t1*sigma_t2,rho_t*sigma_t1*sigma_t3,
                      rho_t*sigma_t2*sigma_t1,sigma_t2^2,rho_t*sigma_t2*sigma_t3,
                      rho_t*sigma_t3*sigma_t1,rho_t*sigma_t3*sigma_t2,sigma_t3^2),byrow=T,3,3)
      cove = matrix(c(sigma_e1^2,rho_e*sigma_e1*sigma_e2,rho_e*sigma_e1*sigma_e3,
                      rho_e*sigma_e2*sigma_e1,sigma_e2^2,rho_e*sigma_e2*sigma_e3,
                      rho_e*sigma_e3*sigma_e1,rho_e*sigma_e3*sigma_e2,sigma_e3^2),byrow=T,3,3)
      naive = matrix(0,N,ncol = 3)
      unive1 = unive2 = unive3 = derived1 = derived2 = derived3 = rep(N,0)
      for(i in 1:N){
        set.seed(123456+i)
        dtaset = gend(n,c,trube,mut,covt,cove,beta_q)
        t = dtaset[,1]; delta = dtaset[,2]
        t1 = dtaset[,3]; t2 = dtaset[,4]; t3 = dtaset[,5]
        x1 = dtaset[,6]; x2 = dtaset[,7]; x3 = dtaset[,8] 
        #naive estimation
        naive_cox = coxph(Surv(t,delta)~x1+x2+x3) 
        naive[i,] = naive_cox$coefficients
        #univariate regression
        univ_cox1 = coxph(Surv(t,delta)~x1)
        unive1[i] = univ_cox1$coefficients
        univ_cox2 = coxph(Surv(t,delta)~x2)
        unive2[i] = univ_cox2$coefficients
        univ_cox3 = coxph(Surv(t,delta)~x3)
        unive3[i] = univ_cox3$coefficients
        #derived estimation
        x2dx1 = x2/x1; x3dx1 = x3/x1
        derive_cox1 = coxph(Surv(t,delta)~x1+x2dx1+x3dx1)
        derived1[i] = derive_cox1$coefficients[1]
        x1dx2 = x1/x2; x3dx2 = x3/x2
        derive_cox2 = coxph(Surv(t,delta)~x2+x1dx2+x3dx2)
        derived2[i] = derive_cox2$coefficients[1]
        x1dx3 = x1/x3; x2dx3 = x2/x3
        derive_cox3 = coxph(Surv(t,delta)~x3+x1dx3+x2dx3)
        derived3[i] = derive_cox3$coefficients[1]
        #print(i)
      }
      #naive
      naive_res = colMeans(naive)
      b1e_n[j,k] = naive_res[1]
      b2e_n[j,k] = naive_res[2]
      b3e_n[j,k] = naive_res[3]
      #univariate
      b1e_un[j,k] = mean(unive1)
      b2e_un[j,k] = mean(unive2)
      b3e_un[j,k] = mean(unive3)
      #derived
      b1e_d[j,k] = mean(derived1)
      b2e_d[j,k] = mean(derived2)
      b3e_d[j,k] = mean(derived3)
      print(c(j,k))
    }
  }
  ##summary: naive
  bias_b1_naive = b1e_n-trube[1]
  bias_b2_naive = b2e_n-trube[2]
  bias_b3_naive = b3e_n-trube[3]
  #univariate
  bias_b1_univ = b1e_un-trube[1]
  bias_b2_univ = b2e_un-trube[2]
  bias_b3_univ = b3e_un-trube[3]
  #derived
  bias_b1_deriv = b1e_d-trube[1]
  bias_b2_deriv = b2e_d-trube[2]
  bias_b3_deriv = b3e_d-trube[3]
  #######
  bias_naive1 = cbind(rho_e_set,bias_b1_naive)
  bias_naive2 = cbind(rho_e_set,bias_b2_naive)
  bias_naive3 = cbind(rho_e_set,bias_b3_naive)
  colnames(bias_naive1) = colnames(bias_naive2) = colnames(bias_naive3) = c("rho_e",paste0("rho_T=",seq(0.1,0.4,by=0.1)))
  rownames(bias_naive1) = rownames(bias_naive2) = rownames(bias_naive3) = rep("Multivariate",4)
  
  bias_univ1 = cbind(rho_e_set,bias_b1_univ)
  bias_univ2 = cbind(rho_e_set,bias_b2_univ)
  bias_univ3 = cbind(rho_e_set,bias_b3_univ)
  colnames(bias_univ1) = colnames(bias_univ2) = colnames(bias_univ3) = c("rho_e",paste0("rho_T=",seq(0.1,0.4,by=0.1)))
  rownames(bias_univ1) = rownames(bias_univ2) = rownames(bias_univ3) = rep("Univariate",4)
  
  
  bias_deriv1 = cbind(rho_e_set,bias_b1_deriv)
  bias_deriv2 = cbind(rho_e_set,bias_b2_deriv)
  bias_deriv3 = cbind(rho_e_set,bias_b3_deriv)
  colnames(bias_deriv1) = colnames(bias_deriv2) = colnames(bias_deriv3) = c("rho_e",paste0("rho_T=",seq(0.1,0.4,by=0.1)))
  rownames(bias_deriv1) = rownames(bias_deriv2) = rownames(bias_deriv3) = rep("Derived",4)
  
  bias_b1_res = rbind(bias_naive1,bias_univ1,bias_deriv1)
  bias_b2_res = rbind(bias_naive2,bias_univ2,bias_deriv2)
  bias_b3_res = rbind(bias_naive3,bias_univ3,bias_deriv3)
  write.csv(round(bias_b1_res,3), paste0('result_cox_bias_b1_',scenario,'.csv'))
  write.csv(round(bias_b2_res,3), paste0('result_cox_bias_b2_',scenario,'.csv'))
  write.csv(round(bias_b3_res,3), paste0('result_cox_bias_b3_',scenario,'.csv'))
}