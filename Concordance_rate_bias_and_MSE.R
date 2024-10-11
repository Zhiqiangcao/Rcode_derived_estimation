###Concordance rate, bias and MSE of naive estimation, univariate estimation, derived estimation
###and RC estimation
###this program can reproduce simulations results in Table 5 of supplementary material
###author: Zhiqiang CAO
###date: 2024.8.15


rm(list=ls())
library(survival) 
library(mvtnorm)
library(MASS)

###############data generation################################
gend = function(n,c,b,mut,covt,cove,beta_q){
  #true variables
  datat = rmvnorm(n=n, mean=mut, sigma=covt)
  u = runif(n)
  t = -log(u)/exp(datat%*%b) #baseline hazard function lambda_0(t)=1
  cen = runif(n, min=0, max=c)
  delta = 1*(t<=cen)    
  time = delta*t+(1-delta)*cen  #observed time
  #error variables
  datae = rmvnorm(n=n, mean=c(0,0,0), sigma=cove) 
  e1 = datae[,1]; e2 = datae[,2]; e3 = datae[,3]; 
  #additive error model
  t1 = datat[,1]; t2 = datat[,2]; t3 = datat[,3];
  #observed variables
  x1 = beta_q[1]*t1+e1; 
  x2 = beta_q[2]*t2+e2
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
scenario = 1 #scenario=1 is for male; scenario=2 is for female

if(scenario == 1){ #X1 = protein, x2 = fat, X3 = energy, male 
  mut = c(9.356683,8.654661,11.937024)
  trube = c(0.172,0.228,0.187)#
  covt=matrix(c(1.05927134,0.050755983,0.33722698,
                0.05075598,0.676466609,0.31167934,
                0.33722698,0.31167934,1.13869850),byrow=T,3,3)
  covt = round(covt,4)
  cove = matrix(c(0.8575957,0.7162641,0.7517087,
                  0.7162641,0.8836357,0.7924088,
                  0.7517087,0.7924088,0.9131195),byrow=T,3,3)
  cove = round(cove,4)
  sigtx = matrix(c(0.38838743,0.021051055,0.093149182,
                   0.01860995,0.280564669,0.086092387,
                   0.12364605,0.129269071,0.314532470),byrow=T,3,3)
  sigtx = round(sigtx,4)
  beta_q = c(0.3667,0.4148,0.2762)
}else if(scenario == 2){ #female
  mut = c(9.245602,7.510393,11.167943)
  trube = c(0.16,0.112,0.074) 
  covt = matrix(c(1.07312719,0.05953603,0.32754194,
                  0.05953603,1.07539124,0.44104340,
                  0.32754194,0.44104340,1.17795383),byrow=T,3,3)
  covt = round(covt,4)
  cove = matrix(c(0.9046743,0.7202257,0.7603176,
                  0.7202257,0.9308771,0.8315039,
                  0.7603176,0.8315039,0.9254783),byrow=T,3,3)
  cove = round(cove,4)
  sigtx = matrix(c(0.31983837,0.015094115,0.082384259,
                   0.01774431,0.272642958,0.110932463,
                   0.09762168,0.111817330,0.296282223),byrow=T,3,3)
  sigtx = round(sigtx,4)
  beta_q = c(0.2980,0.2535,0.2515)
}

if(scenario<2){
  n = 5785  #sample size for male
}else{
  n = 9472  #sample size for female
}

##keep simulation results
naive = naive_p = univ = univ_p = derived = derived_p = rce = rce_p = matrix(0,N,ncol = 3)

#iteration
for(i in 1:N){
  set.seed(10000+i)
  #set.seed(i)
  dtaset = gend(n,c=1.6,trube,mut,covt,cove,beta_q)
  t = dtaset[,1]; delta = dtaset[,2]
  t1 = dtaset[,3]; t2 = dtaset[,4]; t3 = dtaset[,5]
  x1 = dtaset[,6]; x2 = dtaset[,7]; x3 = dtaset[,8]
  #naive estimation
  naive_cox = coxph(Surv(t,delta)~x1+x2+x3)
  naive_summ = summary(naive_cox)$coefficients
  naive[i,] = naive_summ[,1]
  naive_p[i,] = naive_summ[,5]
  #rc estimation
  x = cbind(x1,x2,x3)
  sigxx = cov(x)  #variance-covariance between variables measured with error
  umx = apply(x,2,mean)
  etx = matrix(rep(mut,n),byrow=TRUE,ncol=3)+t(sigtx%*%solve(sigxx)%*%t(x-umx)) 
  m1 = etx[,1]; m2 = etx[,2]; m3 = etx[,3];
  rc_cox = coxph(Surv(t,delta)~m1+m2+m3)
  rc_summ =  summary(rc_cox)$coefficients
  rce[i,] = rc_summ[,1]
  rce_p[i,] = rc_summ[,5]
  #univariate estimation
  univ_cox1 = coxph(Surv(t,delta)~x1)
  univ_summ1 =  summary(univ_cox1)$coefficients
  univ[i,1] = univ_summ1[1,1]
  univ_p[i,1] = univ_summ1[1,5]
  
  univ_cox2 = coxph(Surv(t,delta)~x2)
  univ_summ2 =  summary(univ_cox2)$coefficients
  univ[i,2] = univ_summ2[1,1]
  univ_p[i,2] = univ_summ2[1,5]
  
  univ_cox3 = coxph(Surv(t,delta)~x3)
  univ_summ3 =  summary(univ_cox3)$coefficients
  univ[i,3] = univ_summ3[1,1]
  univ_p[i,3] = univ_summ3[1,5]
  #derived estimation
  x1d = x1/x3; x2d = x2/x3
  derive_cox3 = coxph(Surv(t,delta)~x1d+x2d+x3)
  derive_summ3 =  summary(derive_cox3)$coefficients
  derived[i,3] = derive_summ3[3,1]
  derived_p[i,3] = derive_summ3[3,5]
  
  x1df = x1/x2; x3df = x3/x2
  derive_cox2 = coxph(Surv(t,delta)~x1df+x2+x3df)
  derive_summ2 =  summary(derive_cox2)$coefficients
  derived[i,2] = derive_summ2[2,1]
  derived_p[i,2] = derive_summ2[2,5]
  
  x2dff = x2/x1; x3dff = x3/x1
  derive_cox1 = coxph(Surv(t,delta)~x1+x2dff+x3dff)
  derive_summ1 =  summary(derive_cox1)$coefficients
  derived[i,1] = derive_summ1[1,1]
  derived_p[i,1] = derive_summ1[1,5]
}


###summary estimation
######bias
#naive
naive_res = colMeans(naive)
naive_bias = naive_res-trube
#univariate
univ_res = colMeans(univ)
univ_bias = univ_res-trube
#derived
derive_res = colMeans(derived)
derived_bias = derive_res-trube
#rc
rc_res = colMeans(rce)
rc_bias = rc_res-trube
Bias_res = data.frame(Naive=naive_bias,Univariate=univ_bias,
                      Derived=derived_bias,RC=rc_bias)
row.names(Bias_res) = paste0("beta",1:3)

#####MSE
trubem = matrix(rep(trube,N),byrow=T,ncol=3)
naive_mse = colMeans((naive-trubem)^2)
univ_mse = colMeans((univ-trubem)^2)
derived_mse = colMeans((derived-trubem)^2)
rc_mse = colMeans((rce-trubem)^2)
MSE_res = data.frame(Naive=naive_mse,Univariate=univ_mse,
                     Derived=derived_mse,RC=rc_mse)
row.names(MSE_res) = paste0("beta",1:3)

#####Concordance rate
cal_cr = function(est_p,rce_p){
  #est_p can be naive_p, univ_p and derived_p, which is a N*3 matrix
  cr_res = numeric(3)
  for(j in 1:3){
    #comparing results with RC
    n00 = sum(est_p[,j]>=0.05 & rce_p[,j]>=0.05) #both cannot reject
    n01 = sum(est_p[,j]>=0.05 & rce_p[,j]<0.05)  #naive cannot reject but rc can reject
    n10 = sum(est_p[,j]<0.05 & rce_p[,j]>=0.05)  #naive reject but rc cannot reject
    n11 = sum(est_p[,j]<0.05 & rce_p[,j]<0.05)   #both reject
    cr_res[j] = (n00+n11)/N #ratio of both reject and both cannot reject
  }
  return(cr_res)
}

naive_cr = cal_cr(naive_p,rce_p)
univ_cr = cal_cr(univ_p,rce_p)
derived_cr = cal_cr(derived_p,rce_p)

CR_res = data.frame(Naive=naive_cr,Univariate=univ_cr,Derived=derived_cr)
row.names(CR_res) = paste0("beta",1:3)

#keep results
write.csv(round(CR_res,3), paste0('cr_res_',scenario,'.csv'))
write.csv(round(Bias_res,3), paste0('bias_res_',scenario,'.csv'))
write.csv(round(MSE_res,3), paste0('mse_res_',scenario,'.csv'))