###Bias from bivariate, ridge, univariate and derived estimation based on Cox regression
###this program can reproduce simulation results in Tables 4 of main text
###author: Zhiqiang Cao

rm(list = ls(all = TRUE))
library(mvtnorm)
library(MASS)
library(survival) 
library(glmnet)


###############data generation################################
gend = function(n,c,b,mut,covt,cove){
  b1 = b[1]; b2 = b[2]
  #true variables
  datat = rmvnorm(n=n, mean=mut, sigma=covt)
  t1 = datat[,1]; t2 = datat[,2] 
  u = runif(n)
  t = -log(u)/exp(b1*t1+b2*t2) #baseline hazard function lambda_0(t)=1
  cen = runif(n, min=0, max=c)
  delta = 1*(t<=cen)    
  time = delta*t+(1-delta)*cen  #observed time
  #error variables
  datae = rmvnorm(n=n, mean=c(0,0), sigma=cove) 
  e1 = datae[,1]; e2 = datae[,2]
  #classical error model
  x1 = t1+e1 
  x2 = t2+e2 
  datat = data.frame(
    t = time,
    status = delta,
    t1 = t1,
    t2 = t2,
    x1 = x1,
    x2 = x2
  )
  return(datat)
}


###parameter settings
sigma_x1 = sigma_x2 = 1 #standard error of observed covariates
N = 1000  #simulation times
n = 400   #sample size for each simulation 
rho_1 = 0.5 #reliability index of X1
rho_2 = 0.6 #reliability index of X2
trube = c(1,2)

scenario_set = 1:2
for(scenario in scenario_set){
  if(scenario == 1){
    #range of beta in Table 4 for both male and female are [1.08-2.16]
    mut = c(1,1) #mean vector of true covariates
    c = 0.08     #produce about 50% censoring  
  }else if(scenario == 2){
    mut = c(5,5) #mean vector of true covariates
    c = 5*10^{-7}  #produce about 50% censoring  
  }
  
  sigma_t1 = sigma_x1*rho_1; #square-root of variance of true variable X1
  sigma_t2 = sigma_x2*rho_2; #square-root of variance of true variable X2
  sigma_e1 = sigma_x1*sqrt(1-rho_1^2) #square-root of variance of error variable e1
  sigma_e2 = sigma_x2*sqrt(1-rho_2^2) #square-root of variance of error variable e2
  
  #set of \rho_T
  rho_t_set = c(0.3,0.5,0.7,0.9)
  len_rhot = length(rho_t_set)
  #set of \rho_{\epsilon}
  rho_e_set = seq(0.5,0.9,by=0.1)
  len_rhoe = length(rho_e_set)
  b1e_n = b2e_n = b1e_un = b2e_un = b1e_d = b2e_d = b1e_r = b2e_r = matrix(0,len_rhoe,len_rhot)
  for(k in 1:len_rhot){
    rho_t = rho_t_set[k]
    for(j in 1:len_rhoe){
      rho_e = rho_e_set[j]
      covt = matrix(c(sigma_t1^2,rho_t*sigma_t1*sigma_t2,
                      rho_t*sigma_t1*sigma_t2,sigma_t2^2),2,2)
      cove = matrix(c(sigma_e1^2,rho_e*sigma_e1*sigma_e2,
                      rho_e*sigma_e1*sigma_e2,sigma_e2^2),2,2)
      naive = ridge = matrix(0,N,ncol = 2)
      unive1 = unive2 = derived1 = derived2 = rep(N,0)
      for(i in 1:N){
        set.seed(123456+i)
        dtaset = gend(n,c,trube,mut,covt,cove)
        t = dtaset[,1]; status = dtaset[,2]
        t1 = dtaset[,3]; t2 = dtaset[,4] #true covariates
        x1 = dtaset[,5]; x2 = dtaset[,6] #observed covariates
        #naive estimation
        naive_cox = coxph(Surv(t,status)~x1+x2)
        naive[i,] = naive_cox$coefficients
        #univariate regression
        univ_cox1 = coxph(Surv(t,status)~x1)
        unive1[i] = univ_cox1$coefficients
        univ_cox2 = coxph(Surv(t,status)~x2)
        unive2[i] = univ_cox2$coefficients
        #derived estimation
        x2dx1 = x2/x1; 
        derive_cox1 = coxph(Surv(t,status)~x1+x2dx1)
        derived1[i] = derive_cox1$coefficients[1]
        x1dx2 = x1/x2
        derive_cox2 = coxph(Surv(t,status)~x1dx2+x2)
        derived2[i] = derive_cox2$coefficients[2]
        #ridge regression
        time = t
        y = cbind(time,status)
        x = cbind(x1,x2)
        cvfit = cv.glmnet(x, y, family = "cox", type.measure = "C")
        best_lambda = cvfit$lambda.min
        glmnet_fit = glmnet(x, y, family = "cox", lambda = best_lambda)
        ridge[i,] = coef(glmnet_fit)[1:2]
      }
      naive_res = colMeans(naive)
      b1e_n[j,k] = naive_res[1]
      b2e_n[j,k] = naive_res[2]
      b1e_un[j,k] = mean(unive1)
      b2e_un[j,k] = mean(unive2)
      b1e_d[j,k] = mean(derived1)
      b2e_d[j,k] = mean(derived2)
      ridge_res = colMeans(ridge)
      b1e_r[j,k] = ridge_res[1]
      b2e_r[j,k] = ridge_res[2]
      print(c(j,k))
    }
  }
  ##summary
  bias_b1_naive = b1e_n-trube[1]
  bias_b2_naive = b2e_n-trube[2]
  bias_b1_univ = b1e_un-trube[1]
  bias_b2_univ = b2e_un-trube[2]
  bias_b1_deriv = b1e_d-trube[1]
  bias_b2_deriv = b2e_d-trube[2]
  bias_b1_ridge = b1e_r-trube[1]
  bias_b2_ridge = b2e_r-trube[2]
  #naive
  bias_naive1 = cbind(rho_e_set,bias_b1_naive)
  bias_naive2 = cbind(rho_e_set,bias_b2_naive)
  colnames(bias_naive1) = colnames(bias_naive2) = c("rho_e",paste0("rho_T=",seq(0.3,0.9,by=0.2)))
  rownames(bias_naive1) = rownames(bias_naive2) = rep("Bivariate",5)
  #univariate
  bias_univ1 = cbind(rho_e_set,bias_b1_univ)
  bias_univ2 = cbind(rho_e_set,bias_b2_univ)
  colnames(bias_univ1) = colnames(bias_univ2) = c("rho_e",paste0("rho_T=",seq(0.3,0.9,by=0.2)))
  rownames(bias_univ1) = rownames(bias_univ2) = rep("Univariate",5)
  #derived
  bias_deriv1 = cbind(rho_e_set,bias_b1_deriv)
  bias_deriv2 = cbind(rho_e_set,bias_b2_deriv)
  colnames(bias_deriv1) = colnames(bias_deriv2) = c("rho_e",paste0("rho_T=",seq(0.3,0.9,by=0.2)))
  rownames(bias_deriv1) = rownames(bias_deriv2) = rep("Derived",5)
  #ridge
  bias_ridge1 = cbind(rho_e_set,bias_b1_ridge)
  bias_ridge2 = cbind(rho_e_set,bias_b2_ridge)
  colnames(bias_ridge1) = colnames(bias_ridge2) = c("rho_e",paste0("rho_T=",seq(0.3,0.9,by=0.2)))
  rownames(bias_ridge1) = rownames(bias_ridge2) = rep("ridge",5)
  #summary overall
  bias_b1_res = rbind(bias_naive1,bias_ridge1,bias_univ1,bias_deriv1)
  bias_b2_res = rbind(bias_naive2,bias_ridge2,bias_univ2,bias_deriv2)
  write.csv(round(bias_b1_res,3), paste0('biva_cox_bias_b1_',scenario,'.csv'))
  write.csv(round(bias_b2_res,3), paste0('biva_cox_bias_b2_',scenario,'.csv'))
}