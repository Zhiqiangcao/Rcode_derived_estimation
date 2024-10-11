###Bias from bivariate, univariate and derived three methods of estimation based on linear regression
###this program can reproduce simulations results in Tables 1-2 of main text and Tables S1-S5 of supplementary material
###author: Zhiqiang CAO
###date: 2024.8.6

rm(list = ls(all = TRUE))
library(mvtnorm)
library(MASS)

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
sigma_x1 = sigma_x2 = 1 #variance of observed covariates
N = 1000  #simulation times
n = 400   #sample size for each simulation 
mut = c(6,6) #mean vector of true covariates
#scenario 1-4 is to reproduce simulation results in Table 1 of main text
#scenario 5-8 is to reproduce simulation results in Table 2 of main text
#scenario 9-12 is to reproduce simulation results in Table S5 of supplementary materials
#If we change mut to different values, that is, mut=c(1,1),c(3,3),c(6,6),c(10,10),c(10,20),c(20,10) and c(20,20)
#and focus on scenarios 1-4, then we can reproduce simulation results of Table S1-S4 of supplementary materials
scenario_set = 1:4
for(scenario in scenario_set){
  if(scenario == 1){
    trube = c(1,3)
    rho_1 = 0.5
    rho_2 = 0.9
  }else if(scenario == 2){
    trube = c(1,3)
    rho_1 = 0.5
    rho_2 = 0.6
  }else if(scenario == 3){
    trube = c(3,1)
    rho_1 = 0.5
    rho_2 = 0.9
  }else if(scenario == 4){
    trube = c(3,1)
    rho_1 = 0.5
    rho_2 = 0.6
  }else if(scenario == 5){
    trube = c(1,-3)
    rho_1 = 0.5
    rho_2 = 0.9
  }else if(scenario == 6){
    trube = c(1,-3)
    rho_1 = 0.5
    rho_2 = 0.6
  }else if(scenario == 7){
    trube = c(3,-1)
    rho_1 = 0.5
    rho_2 = 0.9
  }else if(scenario == 8){
    trube = c(3,-1)
    rho_1 = 0.5
    rho_2 = 0.6
  }else if(scenario == 9){
    trube = c(1,1.5)
    rho_1 = 0.5
    rho_2 = 0.9
  }else if(scenario == 10){
    trube = c(1,1.5)
    rho_1 = 0.5
    rho_2 = 0.6
  }else if(scenario == 11){
    trube = c(1.5,1)
    rho_1 = 0.5
    rho_2 = 0.9
  }else if(scenario == 12){
    trube = c(1.5,1)
    rho_1 = 0.5
    rho_2 = 0.6
  }
  
  sigma_t1 = rho_1; #square-root of variance of true variable X1
  sigma_t2 = rho_2; #square-root of variance of true variable X2
  sigma_e1 = sqrt(sigma_x1^2-rho_1^2) #square-root of variance of error variable e1
  sigma_e2 = sqrt(sigma_x2^2-rho_2^2) #square-root of variance of error variable e2
  
  #set of \rho_T
  rho_t_set = c(0.3,0.5,0.7,0.9)
  len_rhot = length(rho_t_set)
  #set of \rho_{\epsilon}
  rho_e_set = seq(0.5,0.9,by=0.1)
  len_rhoe = length(rho_e_set)
  
  b1e_n = b2e_n = b1e_un = b2e_un = b1e_d = b2e_d = matrix(0,len_rhoe,len_rhot)
  for(k in 1:len_rhot){
    rho_t = rho_t_set[k]
    for(j in 1:len_rhoe){
      rho_e = rho_e_set[j]
      covt = matrix(c(sigma_t1^2,rho_t*sigma_t1*sigma_t2,
                      rho_t*sigma_t1*sigma_t2,sigma_t2^2),2,2)
      cove = matrix(c(sigma_e1^2,rho_e*sigma_e1*sigma_e2,
                      rho_e*sigma_e1*sigma_e2,sigma_e2^2),2,2)
      naive = matrix(0,N,ncol = 2)
      unive1 = unive2 = derived1 = derived2 = rep(N,0)
      for(i in 1:N){
        set.seed(123456+i)
        dtaset = gend(n,trube,mut,covt,cove)
        y = dtaset[,1]; 
        t1 = dtaset[,2]; t2 = dtaset[,3] #true covariates
        x1 = dtaset[,4]; x2 = dtaset[,5] #observed covariates
        #naive estimation
        naive_lm = lm(y~x1+x2)
        naive[i,] = naive_lm$coefficients[2:3]
        #univariate regression
        univ_lm1 = lm(y~x1)
        unive1[i] = univ_lm1$coefficients[2]
        univ_lm2 = lm(y~x2)
        unive2[i] = univ_lm2$coefficients[2]
        #derived estimation
        x2dx1 = x2/x1
        derive_lm1 = lm(y~x1+x2dx1)
        derived1[i] = derive_lm1$coefficients[2]
        x1dx2 = x1/x2
        derive_lm2 = lm(y~x1dx2+x2)
        derived2[i] = derive_lm2$coefficients[3]
        #print(i)
      }
      naive_res = colMeans(naive)
      b1e_n[j,k] = naive_res[1]
      b2e_n[j,k] = naive_res[2]
      b1e_un[j,k] = mean(unive1)
      b2e_un[j,k] = mean(unive2)
      b1e_d[j,k] = mean(derived1)
      b2e_d[j,k] = mean(derived2)
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
  bias_naive1 = cbind(rho_e_set,bias_b1_naive)
  bias_naive2 = cbind(rho_e_set,bias_b2_naive)
  colnames(bias_naive1) = colnames(bias_naive2) = c("rho_e",paste0("rho_T=",seq(0.3,0.9,by=0.2)))
  rownames(bias_naive1) = rownames(bias_naive2) = rep("Bivariate",5)
  bias_univ1 = cbind(rho_e_set,bias_b1_univ)
  bias_univ2 = cbind(rho_e_set,bias_b2_univ)
  colnames(bias_univ1) = colnames(bias_univ2) = c("rho_e",paste0("rho_T=",seq(0.3,0.9,by=0.2)))
  rownames(bias_univ1) = rownames(bias_univ2) = rep("Univariate",5)
  bias_deriv1 = cbind(rho_e_set,bias_b1_deriv)
  bias_deriv2 = cbind(rho_e_set,bias_b2_deriv)
  colnames(bias_deriv1) = colnames(bias_deriv2) = c("rho_e",paste0("rho_T=",seq(0.3,0.9,by=0.2)))
  rownames(bias_deriv1) = rownames(bias_deriv2) = rep("Derived",5)
  
  bias_b1_res = rbind(bias_naive1,bias_univ1,bias_deriv1)
  bias_b2_res = rbind(bias_naive2,bias_univ2,bias_deriv2)
  write.csv(round(bias_b1_res,3), paste0('result_lm_bias_b1_',scenario,'.csv'))
  write.csv(round(bias_b2_res,3), paste0('result_lm_bias_b2_',scenario,'.csv'))
}