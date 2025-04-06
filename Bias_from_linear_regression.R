###Bias from bivariate, ridge, univariate and derived three methods of estimation based on linear regression
###this program can reproduce simulations results in Tables 2-3 of main text and Tables S1-S2 of supplementary material
###author: Zhiqiang Cao

rm(list = ls(all = TRUE))
library(mvtnorm)
library(MASS)
library(glmnet)

###############data generation################################
gend = function(n,b,mut,covt,cove){
  #true regression coefficients
  b1 = b[1]; b2 = b[2]
  #true variables
  datat = rmvnorm(n=n, mean=mut, sigma=covt)
  t1 = datat[,1]; t2 = datat[,2] 
  y = b1*t1+b2*t2+rnorm(n,0,sd=0.1) #plus a error term
  #error variables
  datae = rmvnorm(n=n, mean=c(0,0), sigma=cove) 
  e1 = datae[,1]; e2 = datae[,2]
  #classical error model
  x1 = t1+e1 
  x2 = t2+e2 
  datat = data.frame(
    y = y,
    t1 = t1,
    t2 = t2,
    x1 = x1,
    x2 = x2
  )
  return(datat)
}


###parameter settings
sigma_x1 = sigma_x2 = 1 #standard error of observed covariates
#set sigma_x1=sigma_x2 = sqrt(10) or 10 can reproduce simulation results in Tables R2 and R3 of response letter 
N = 1000  #simulation times
n = 400   #sample size for each simulation 
rho_1 = 0.5 #reliability index of X1
rho_2 = 0.6 #reliability index of X2
trube = c(1,2) #this settings can reproduce simulation results in Table 2 of main text
#change trube = c(2,1) can reproduce simulation results in Table 3 of main text
#change trube = c(1,-2) can reproduce simulation results in Table S1 of Web Appendix E
#change trube = c(1,1.5) can reproduce simulation results in Table S2 of Web Appendix E
scenario_set = 1:5
for(scenario in scenario_set){
  if(scenario == 1){
    mut = c(1,1) 
  }else if(scenario == 2){
    mut = c(5,5) 
  }else if(scenario == 3){
    mut = c(5,10)
  }else if(scenario == 4){
    mut = c(10,5)
  }else if(scenario == 5){
    mut = c(10,10)
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
        dtaset = gend(n,trube,mut,covt,cove)
        y = dtaset[,1]; 
        t1 = dtaset[,2]; t2 = dtaset[,3] #true covariates
        x1 = dtaset[,4]; x2 = dtaset[,5] #observed covariates
        x1s = scale(x1); x2s = scale(x2)
        #naive estimation
        naive_lm = lm(y~x1+x2)
        naive[i,] = naive_lm$coefficients[2:3]
        #univariate regression
        univ_lm1 = lm(y~x1)
        unive1[i] = univ_lm1$coefficients[2]
        univ_lm2 = lm(y~x2)
        unive2[i] = univ_lm2$coefficients[2]
        #derived estimation
        x2dx1 = x2/x1; x2dx1s = scale(x2dx1)
        derive_lm1 = lm(y~x1+x2dx1)
        derived1[i] = derive_lm1$coefficients[2]
        x1dx2 = x1/x2; x1dx2s = scale(x1dx2)
        derive_lm2 = lm(y~x1dx2+x2)
        derived2[i] = derive_lm2$coefficients[3]
        #ridge regression
        x = cbind(x1,x2)
        cv_model = cv.glmnet(x, y, alpha = 0)
        best_lambda = cv_model$lambda.min
        best_model = glmnet(x, y, alpha = 0, lambda = best_lambda, standardize = FALSE)
        ridge[i,] = coef(best_model)[2:3]
        #print(i)
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
  ##summary simulation results
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
  #ridge
  bias_ridge1 = cbind(rho_e_set,bias_b1_ridge)
  bias_ridge2 = cbind(rho_e_set,bias_b2_ridge)
  colnames(bias_ridge1) = colnames(bias_ridge2) = c("rho_e",paste0("rho_T=",seq(0.3,0.9,by=0.2)))
  rownames(bias_ridge1) = rownames(bias_ridge2) = rep("ridge",5)
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
  #overall rbind
  bias_b1_res = rbind(bias_naive1,bias_ridge1,bias_univ1,bias_deriv1)
  bias_b2_res = rbind(bias_naive2,bias_ridge2,bias_univ2,bias_deriv2)
  write.csv(round(bias_b1_res,3), paste0('biva_lm_bias_b1_',scenario,'.csv'))
  write.csv(round(bias_b2_res,3), paste0('biva_lm_bias_b2_',scenario,'.csv'))
}