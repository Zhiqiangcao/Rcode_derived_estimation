# functions used in the simulation studies


### data generation for bivariate linear models 
# parameter description 
# n: sample size
# b: regression parameter 
# mut: mean of true covariates 
# covt: variance-covariance matrix of true covariates 
# cove: variance-covariance matrix of error covariates
gend_biva_linear = function(n, b, mut, covt, cove) {
  b1 = b[1]
  b2 = b[2]
  # true variables
  datat = rmvnorm(n = n, mean = mut, sigma = covt)
  t1 = datat[, 1]
  t2 = datat[, 2]
  y = b1 * t1 + b2 * t2 + rnorm(n, 0, sd = 0.1) 
  # error variables
  datae = rmvnorm(n = n, mean = c(0, 0), sigma = cove)
  e1 = datae[, 1]
  e2 = datae[, 2]
  # classical error model
  x1 = t1 + e1
  x2 = t2 + e2
  datat = data.frame(y = y, t1 = t1, t2 = t2, x1 = x1, x2 = x2)
  return(datat)
}

### data generation for bivariate Cox model 
# parameter description 
# n: sample size
# c: parameter used to generate censoring variables
# b: regression parameter 
# mut: mean of true covariates 
# covt: variance-covariance matrix of true covariates 
# cove: variance-covariance matrix of error covariates
gend_biva_cox = function(n, c, b, mut, covt, cove) {
  b1 = b[1]
  b2 = b[2]
  # true variables
  datat = rmvnorm(n = n, mean = mut, sigma = covt)
  t1 = datat[, 1]
  t2 = datat[, 2]
  u = runif(n)
  t = -log(u)/exp(b1 * t1 + b2 * t2)  #baseline hazard function lambda_0(t)=1
  cen = runif(n, min = 0, max = c)
  delta = 1 * (t <= cen)
  time = delta * t + (1 - delta) * cen  
  # error variables
  datae = rmvnorm(n = n, mean = c(0, 0), sigma = cove)
  e1 = datae[, 1]
  e2 = datae[, 2]
  # classical error model
  x1 = t1 + e1
  x2 = t2 + e2
  datat = data.frame(t = time, status = delta, t1 = t1, t2 = t2, x1 = x1, x2 = x2)
  return(datat)
}

### data generation for multivariate linear model 
# parameter description 
# n: sample size
# b: regression parameter 
# mut: mean of true covariates 
# covt: variance-covariance matrix of true covariates 
# cove: variance-covariance matrix of error covariates
# beta_q: bias intercept parameter in additive error model
gend_mult_linear = function(n, b, mut, covt, cove, beta_q) {
  b1 = b[1]
  b2 = b[2]
  b3 = b[3]
  # true variables
  datat = rmvnorm(n = n, mean = mut, sigma = covt)
  t1 = datat[, 1]
  t2 = datat[, 2]
  t3 = datat[, 3]
  y = b1 * t1 + b2 * t2 + b3 * t3 + rnorm(n, 0, sd = 0.1)  
  # error variables
  datae = rmvnorm(n = n, mean = c(0, 0, 0), sigma = cove)
  e1 = datae[, 1]
  e2 = datae[, 2]
  e3 = datae[, 3]
  # additive error model
  x1 = beta_q[1] * t1 + e1
  x2 = beta_q[2] * t2 + e2
  x3 = beta_q[3] * t3 + e3
  datat = data.frame(y = y, t1 = t1, t2 = t2, t3 = t3, x1 = x1, x2 = x2, x3 = x3)
  return(datat)
}

### data generation for multivariate Cox model 
# parameter description 
# n: sample size
# c: parameter used to generate censoring variables
# b: regression parameter 
# mut: mean of true covariates 
# covt: variance-covariance matrix of true covariates 
# cove: variance-covariance matrix of error covariates
# beta_q: bias intercept parameter in additive error model
gend_mult_cox = function(n, c, b, mut, covt, cove, beta_q) {
  b1 = b[1]
  b2 = b[2]
  b3 = b[3]
  # true variables
  datat = rmvnorm(n = n, mean = mut, sigma = covt)
  t1 = datat[, 1]
  t2 = datat[, 2]
  t3 = datat[, 3]
  u = runif(n)
  t = -log(u)/exp(b1 * t1 + b2 * t2 + b3 * t3)  #baseline hazard function lambda_0(t)=1
  cen = runif(n, min = 0, max = c)
  delta = 1 * (t <= cen)
  time = delta * t + (1 - delta) * cen 
  # error variables
  datae = rmvnorm(n = n, mean = c(0, 0, 0), sigma = cove)
  e1 = datae[, 1]
  e2 = datae[, 2]
  e3 = datae[, 3]
  # additive error model
  x1 = beta_q[1] * t1 + e1
  x2 = beta_q[2] * t2 + e2
  x3 = beta_q[3] * t3 + e3
  datat = data.frame(t = time, status = delta, t1 = t1, t2 = t2, t3 = t3, x1 = x1, x2 = x2,
                     x3 = x3)
  return(datat)
}

# generate simulated data based on a multiplicative error model
# parameter description 
# n: sample size
# b: regression parameter 
# mut: mean of true covariates 
# covt: variance-covariance matrix of true covariates 
# cove: variance-covariance matrix of error covariates
gend_multiplicative_error_data = function(n, b, mut, covt, cove) {
  b1 = b[1]
  b2 = b[2]
  # true variables
  datat = rmvnorm(n = n, mean = mut, sigma = covt)
  t1 = datat[, 1]
  t2 = datat[, 2]
  y = b1 * t1 + b2 * t2 + rnorm(n, 0, sd = 0.1)  
  # error variables
  datae = rmvnorm(n = n, mean = c(0, 0), sigma = cove)
  e1 = datae[, 1]
  e2 = datae[, 2]
  # multiplicative error model
  x1 = t1 * abs(e1)
  x2 = t2 * abs(e2)
  datat = data.frame(y = y, t1 = t1, t2 = t2, x1 = x1, x2 = x2)
  return(datat)
}

### bivariate linear regression, which can reproduce results in Tables 2-3 of
### manuscript and Tables S1-S3, S8 of Web Appendix 
# parameter description
# trube: true regression coefficients 
# sigma_x1: standard error of observed covariate x1 
# sigma_x2: standard error of observed covariate x2 
# error_type: error_type=1 means classical error model; error_type=2 means multiplicative error model 
# ridge_reg: including ridge regression or not
est_biva_linear <- function(trube, sigma_x1, sigma_x2, error_type = 1, ridge_reg = TRUE) {
  if (error_type == 1)
    scenario_set = 1:5 else scenario_set = 1:2
    b1_biva = b2_biva = b1_ridge = b2_ridge = b1_uni = b2_uni = b1_derive = b2_derive = NULL
    for (scenario in scenario_set) {
      if (scenario == 1) {
        mut = c(1, 1)
      } else if (scenario == 2) {
        mut = c(5, 5)
      } else if (scenario == 3) {
        mut = c(5, 10)
      } else if (scenario == 4) {
        mut = c(10, 5)
      } else if (scenario == 5) {
        mut = c(10, 10)
      }
      
      sigma_t1 = sigma_x1 * rho_1  #square-root of variance of true variable X1
      sigma_t2 = sigma_x2 * rho_2  #square-root of variance of true variable X2
      sigma_e1 = sigma_x1 * sqrt(1 - rho_1^2)  #square-root of variance of error variable e1
      sigma_e2 = sigma_x2 * sqrt(1 - rho_2^2)  #square-root of variance of error variable e2
      
      # set of \rho_T
      rho_t_set = c(0.3, 0.5, 0.7, 0.9)
      len_rhot = length(rho_t_set)
      # set of \rho_{\epsilon}
      rho_e_set = seq(0.5, 0.9, by = 0.1)
      len_rhoe = length(rho_e_set)
      
      b1e_n = b2e_n = b1e_un = b2e_un = b1e_d = b2e_d = b1e_r = b2e_r = matrix(0, len_rhoe,
                                                                               len_rhot)
      for (k in 1:len_rhot) {
        rho_t = rho_t_set[k]
        for (j in 1:len_rhoe) {
          rho_e = rho_e_set[j]
          covt = matrix(c(sigma_t1^2, rho_t * sigma_t1 * sigma_t2, rho_t * sigma_t1 *
                            sigma_t2, sigma_t2^2), 2, 2)
          cove = matrix(c(sigma_e1^2, rho_e * sigma_e1 * sigma_e2, rho_e * sigma_e1 *
                            sigma_e2, sigma_e2^2), 2, 2)
          naive = ridge = matrix(0, N, ncol = 2)
          unive1 = unive2 = derived1 = derived2 = rep(N, 0)
          for (i in 1:N) {
            set.seed(123456 + i)
            if (error_type == 1) { #classical error model
              dtaset = gend_biva_linear(n, trube, mut, covt, cove)
            } else { #multiplicative error model 
              dtaset = gend_multiplicative_error_data(n, trube, mut, covt, cove)
            }
            y = dtaset[, 1]
            t1 = dtaset[, 2]
            t2 = dtaset[, 3]  #true covariates
            x1 = dtaset[, 4]
            x2 = dtaset[, 5]  #observed covariates
            # naive estimation
            naive_lm = lm(y ~ x1 + x2)
            naive[i, ] = naive_lm$coefficients[2:3]
            # univariate regression
            univ_lm1 = lm(y ~ x1)
            unive1[i] = univ_lm1$coefficients[2]
            univ_lm2 = lm(y ~ x2)
            unive2[i] = univ_lm2$coefficients[2]
            # derived estimation
            x2dx1 = x2/x1
            derive_lm1 = lm(y ~ x1 + x2dx1)
            derived1[i] = derive_lm1$coefficients[2]
            x1dx2 = x1/x2
            derive_lm2 = lm(y ~ x1dx2 + x2)
            derived2[i] = derive_lm2$coefficients[3]
            # including ridge regression or not
            if (ridge_reg == TRUE) {
              x = cbind(x1, x2)
              cv_model = cv.glmnet(x, y, alpha = 0)
              best_lambda = cv_model$lambda.min
              best_model = glmnet(x, y, alpha = 0, lambda = best_lambda, standardize = FALSE)
              ridge[i, ] = coef(best_model)[2:3]
            }
          }
          naive_res = colMeans(naive)
          b1e_n[j, k] = naive_res[1]
          b2e_n[j, k] = naive_res[2]
          b1e_un[j, k] = mean(unive1)
          b2e_un[j, k] = mean(unive2)
          b1e_d[j, k] = mean(derived1)
          b2e_d[j, k] = mean(derived2)
          ridge_res = colMeans(ridge)
          b1e_r[j, k] = ridge_res[1]
          b2e_r[j, k] = ridge_res[2]
          print(c(j, k))
        }
      }
      ## summary simulation results
      bias_b1_naive = b1e_n - trube[1]
      bias_b2_naive = b2e_n - trube[2]
      bias_b1_univ = b1e_un - trube[1]
      bias_b2_univ = b2e_un - trube[2]
      bias_b1_deriv = b1e_d - trube[1]
      bias_b2_deriv = b2e_d - trube[2]
      bias_b1_ridge = b1e_r - trube[1]
      bias_b2_ridge = b2e_r - trube[2]
      # naive
      bias_naive1 = cbind(rho_e_set, bias_b1_naive)
      bias_naive2 = cbind(rho_e_set, bias_b2_naive)
      colnames(bias_naive1) = colnames(bias_naive2) = c("rho_e", paste0("rho_T=", seq(0.3,
                                                                                      0.9, by = 0.2)))
      rownames(bias_naive1) = rownames(bias_naive2) = rep("Bivariate", 5)
      b1_biva = rbind(b1_biva, bias_naive1)
      b2_biva = rbind(b2_biva, bias_naive2)
      # ridge
      bias_ridge1 = cbind(rho_e_set, bias_b1_ridge)
      bias_ridge2 = cbind(rho_e_set, bias_b2_ridge)
      colnames(bias_ridge1) = colnames(bias_ridge2) = c("rho_e", paste0("rho_T=", seq(0.3,
                                                                                      0.9, by = 0.2)))
      rownames(bias_ridge1) = rownames(bias_ridge2) = rep("ridge", 5)
      b1_ridge = rbind(b1_ridge, bias_ridge1)
      b2_ridge = rbind(b2_ridge, bias_ridge2)
      # univariate
      bias_univ1 = cbind(rho_e_set, bias_b1_univ)
      bias_univ2 = cbind(rho_e_set, bias_b2_univ)
      colnames(bias_univ1) = colnames(bias_univ2) = c("rho_e", paste0("rho_T=", seq(0.3,
                                                                                    0.9, by = 0.2)))
      rownames(bias_univ1) = rownames(bias_univ2) = rep("Univariate", 5)
      b1_uni = rbind(b1_uni, bias_univ1)
      b2_uni = rbind(b2_uni, bias_univ2)
      # derived
      bias_deriv1 = cbind(rho_e_set, bias_b1_deriv)
      bias_deriv2 = cbind(rho_e_set, bias_b2_deriv)
      colnames(bias_deriv1) = colnames(bias_deriv2) = c("rho_e", paste0("rho_T=", seq(0.3,
                                                                                      0.9, by = 0.2)))
      rownames(bias_deriv1) = rownames(bias_deriv2) = rep("Derived", 5)
      b1_derive = rbind(b1_derive, bias_deriv1)
      b2_derive = rbind(b2_derive, bias_deriv2)
    }
    # summary
    b1_Univariate = b1_uni[3, ]
    b2_Univariate = b2_uni[3, ]
    names(b1_Univariate) = names(b2_Univariate) = "Univariate"
    if (error_type == 1) {
      if (ridge_reg == TRUE) {
        b1_final_res = rbind(b1_biva[1:5, ], b1_ridge[1:5, ], b1_Univariate, b1_derive)
        b2_final_res = rbind(b2_biva[1:5, ], b2_ridge[1:5, ], b2_Univariate, b2_derive)
      } else {
        b1_final_res = rbind(b1_biva[1:5, ], b1_Univariate, b1_derive)
        b2_final_res = rbind(b2_biva[1:5, ], b2_Univariate, b2_derive)
      }
    } else {
      b1_final_res = rbind(b1_biva[1:5, ], b1_uni[3, ], b1_derive[1:5, ], 
      b1_biva[6:10, ], b1_uni[8, ], b1_derive[6:10, ])
      b2_final_res = rbind(b2_biva[1:5, ], b2_uni[3, ], b2_derive[1:5, ], 
      b2_biva[6:10, ], b2_uni[8, ], b2_derive[6:10, ])
    }
    res = cbind(b1_final_res, b2_final_res[, -1])
    return(res)
}

# bivariate Cox regression, which can reproduce results in Table 4 of manuscript
est_biva_cox <- function(trube, sigma_x1, sigma_x2) {
  b1_biva = b2_biva = b1_ridge = b2_ridge = b1_uni = b2_uni = b1_derive = b2_derive = NULL
  scenario_set = 1:2
  for (scenario in scenario_set) {
    if (scenario == 1) {
      # range of beta in Table 4 for both male and female are [1.08-2.16]
      mut = c(1, 1)
      c = 0.08  #produce about 50% censoring  
    } else if (scenario == 2) {
      mut = c(5, 5)
      c = 5 * 10^{-7}  #produce about 50% censoring  
    }
    sigma_t1 = sigma_x1 * rho_1  #square-root of variance of true variable X1
    sigma_t2 = sigma_x2 * rho_2  #square-root of variance of true variable X2
    sigma_e1 = sigma_x1 * sqrt(1 - rho_1^2)  #square-root of variance of error variable e1
    sigma_e2 = sigma_x2 * sqrt(1 - rho_2^2)  #square-root of variance of error variable e2
    
    # set of \rho_T
    rho_t_set = c(0.3, 0.5, 0.7, 0.9)
    len_rhot = length(rho_t_set)
    # set of \rho_{\epsilon}
    rho_e_set = seq(0.5, 0.9, by = 0.1)
    len_rhoe = length(rho_e_set)
    b1e_n = b2e_n = b1e_un = b2e_un = b1e_d = b2e_d = b1e_r = b2e_r = matrix(0, len_rhoe,
                                                                             len_rhot)
    for (k in 1:len_rhot) {
      rho_t = rho_t_set[k]
      for (j in 1:len_rhoe) {
        rho_e = rho_e_set[j]
        covt = matrix(c(sigma_t1^2, rho_t * sigma_t1 * sigma_t2, rho_t * sigma_t1 *
                          sigma_t2, sigma_t2^2), 2, 2)
        cove = matrix(c(sigma_e1^2, rho_e * sigma_e1 * sigma_e2, rho_e * sigma_e1 *
                          sigma_e2, sigma_e2^2), 2, 2)
        naive = ridge = matrix(0, N, ncol = 2)
        unive1 = unive2 = derived1 = derived2 = rep(N, 0)
        for (i in 1:N) {
          set.seed(123456 + i)
          dtaset = gend_biva_cox(n, c, trube, mut, covt, cove)
          t = dtaset[, 1]
          status = dtaset[, 2]
          t1 = dtaset[, 3]
          t2 = dtaset[, 4]  #true covariates
          x1 = dtaset[, 5]
          x2 = dtaset[, 6]  #observed covariates
          # naive estimation
          naive_cox = coxph(Surv(t, status) ~ x1 + x2)
          naive[i, ] = naive_cox$coefficients
          # univariate regression
          univ_cox1 = coxph(Surv(t, status) ~ x1)
          unive1[i] = univ_cox1$coefficients
          univ_cox2 = coxph(Surv(t, status) ~ x2)
          unive2[i] = univ_cox2$coefficients
          # derived estimation
          x2dx1 = x2/x1
          derive_cox1 = coxph(Surv(t, status) ~ x1 + x2dx1)
          derived1[i] = derive_cox1$coefficients[1]
          x1dx2 = x1/x2
          derive_cox2 = coxph(Surv(t, status) ~ x1dx2 + x2)
          derived2[i] = derive_cox2$coefficients[2]
          # ridge regression
          time = t
          y = cbind(time, status)
          x = cbind(x1, x2)
          cvfit = cv.glmnet(x, y, family = "cox", type.measure = "C")
          best_lambda = cvfit$lambda.min
          glmnet_fit = glmnet(x, y, family = "cox", lambda = best_lambda)
          ridge[i, ] = coef(glmnet_fit)[1:2]
        }
        naive_res = colMeans(naive)
        b1e_n[j, k] = naive_res[1]
        b2e_n[j, k] = naive_res[2]
        b1e_un[j, k] = mean(unive1)
        b2e_un[j, k] = mean(unive2)
        b1e_d[j, k] = mean(derived1)
        b2e_d[j, k] = mean(derived2)
        ridge_res = colMeans(ridge)
        b1e_r[j, k] = ridge_res[1]
        b2e_r[j, k] = ridge_res[2]
        print(c(j, k))
      }
    }
    ## summary
    bias_b1_naive = b1e_n - trube[1]
    bias_b2_naive = b2e_n - trube[2]
    bias_b1_univ = b1e_un - trube[1]
    bias_b2_univ = b2e_un - trube[2]
    bias_b1_deriv = b1e_d - trube[1]
    bias_b2_deriv = b2e_d - trube[2]
    bias_b1_ridge = b1e_r - trube[1]
    bias_b2_ridge = b2e_r - trube[2]
    # naive
    bias_naive1 = cbind(rho_e_set, bias_b1_naive)
    bias_naive2 = cbind(rho_e_set, bias_b2_naive)
    colnames(bias_naive1) = colnames(bias_naive2) = c("rho_e", paste0("rho_T=", seq(0.3,
                                                                                    0.9, by = 0.2)))
    rownames(bias_naive1) = rownames(bias_naive2) = rep("Bivariate", 5)
    b1_biva = rbind(b1_biva, bias_naive1)
    b2_biva = rbind(b2_biva, bias_naive2)
    # ridge
    bias_ridge1 = cbind(rho_e_set, bias_b1_ridge)
    bias_ridge2 = cbind(rho_e_set, bias_b2_ridge)
    colnames(bias_ridge1) = colnames(bias_ridge2) = c("rho_e", paste0("rho_T=", seq(0.3,
                                                                                    0.9, by = 0.2)))
    rownames(bias_ridge1) = rownames(bias_ridge2) = rep("ridge", 5)
    b1_ridge = rbind(b1_ridge, bias_ridge1)
    b2_ridge = rbind(b2_ridge, bias_ridge2)
    # univariate
    bias_univ1 = cbind(rho_e_set, bias_b1_univ)
    bias_univ2 = cbind(rho_e_set, bias_b2_univ)
    colnames(bias_univ1) = colnames(bias_univ2) = c("rho_e", paste0("rho_T=", seq(0.3,
                                                                                  0.9, by = 0.2)))
    rownames(bias_univ1) = rownames(bias_univ2) = rep("Univariate", 5)
    b1_uni = rbind(b1_uni, bias_univ1)
    b2_uni = rbind(b2_uni, bias_univ2)
    # derived
    bias_deriv1 = cbind(rho_e_set, bias_b1_deriv)
    bias_deriv2 = cbind(rho_e_set, bias_b2_deriv)
    colnames(bias_deriv1) = colnames(bias_deriv2) = c("rho_e", paste0("rho_T=", seq(0.3,
                                                                                    0.9, by = 0.2)))
    rownames(bias_deriv1) = rownames(bias_deriv2) = rep("Derived", 5)
    b1_derive = rbind(b1_derive, bias_deriv1)
    b2_derive = rbind(b2_derive, bias_deriv2)
  }
  # univariate estimation results for each \rho_{\epsilon} are almost the same
  b1_final_res = rbind(b1_biva[1:5, ], b1_ridge[1:5, ], b1_uni[3, ], b1_derive[1:5, ],
                       b1_biva[6:10, ], b1_ridge[6:10, ], b1_uni[8, ], b1_derive[6:10, ])
  b2_final_res = rbind(b2_biva[1:5, ], b2_ridge[1:5, ], b2_uni[3, ], b2_derive[1:5, ],
                       b2_biva[6:10, ], b2_ridge[6:10, ], b2_uni[8, ], b2_derive[6:10, ])
  res = cbind(b1_final_res, b2_final_res[, -1])
  return(res)
}

# multivariate linear or Cox regression, this program can reproduce results in Tables
# S4-S5 of Web Appendix
est_mult_regression <- function(linear = TRUE) {
  scenario_set = 1:2
  b1_biva = b2_biva = b3_biva = b1_uni = b2_uni = b3_uni = b1_derive = b2_derive = b3_derive = NULL
  for (scenario in scenario_set) {
    if (scenario == 1) {
      # for male
      trube = c(0.172, 0.228, 0.187)  #c(1,1.33,1.09)
      mut = c(9.356683, 8.654661, 11.937024)
      c = 0.005  #to produce 50% censoring rate
      beta_q = c(0.3667, 0.4148, 0.2762)
      sigma_t1 = sqrt(1.05927134)
      sigma_t2 = sqrt(0.676466609)
      sigma_t3 = sqrt(1.1386985)
      sigma_e1 = sqrt(0.8575957)
      sigma_e2 = sqrt(0.8836357)
      sigma_e3 = sqrt(0.9131195)
    } else if (scenario == 2) {
      # for female
      trube = c(0.16, 0.112, 0.074)  #c(2.16,1,1.51)
      mux = c(9.245602, 7.510393, 11.167943)
      c = 0.068  #to produce about 50% censoring rate
      beta_q = c(0.298, 0.2535, 0.2515)
      sigma_t1 = sqrt(1.07312719)
      sigma_t2 = sqrt(1.07539124)
      sigma_t3 = sqrt(1.17795383)
      sigma_e1 = sqrt(0.9046743)
      sigma_e2 = sqrt(0.9308771)
      sigma_e3 = sqrt(0.9254783)
    }
    
    rho_t_set = seq(0.1, 0.4, by = 0.1)
    len_rhot = length(rho_t_set)
    rho_e_set = c(0.7, 0.8, 0.9, 0.95)
    len_rhoe = length(rho_e_set)
    b1e_n = b2e_n = b3e_n = b1e_un = b2e_un = b3e_un = b1e_d = b2e_d = b3e_d = matrix(0,
                                                                                      len_rhoe, len_rhot)
    for (k in 1:len_rhot) {
      rho_t = rho_t_set[k]
      for (j in 1:len_rhoe) {
        rho_e = rho_e_set[j]
        covt = matrix(c(sigma_t1^2, rho_t * sigma_t1 * sigma_t2, rho_t * sigma_t1 *
                          sigma_t3, rho_t * sigma_t2 * sigma_t1, sigma_t2^2, rho_t * sigma_t2 * sigma_t3,
                        rho_t * sigma_t3 * sigma_t1, rho_t * sigma_t3 * sigma_t2, sigma_t3^2),
                      byrow = T, 3, 3)
        cove = matrix(c(sigma_e1^2, rho_e * sigma_e1 * sigma_e2, rho_e * sigma_e1 *
                          sigma_e3, rho_e * sigma_e2 * sigma_e1, sigma_e2^2, rho_e * sigma_e2 * sigma_e3,
                        rho_e * sigma_e3 * sigma_e1, rho_e * sigma_e3 * sigma_e2, sigma_e3^2),
                      byrow = T, 3, 3)
        naive = matrix(0, N, ncol = 3)
        unive1 = unive2 = unive3 = derived1 = derived2 = derived3 = rep(N, 0)
        for (i in 1:N) {
          set.seed(123456 + i)
          if (linear == TRUE) {
            # multivariate linear regression
            dtaset = gend_mult_linear(n, trube, mut, covt, cove, beta_q)
            y = dtaset[, 1]
            t1 = dtaset[, 2]
            t2 = dtaset[, 3]
            t3 = dtaset[, 4]
            x1 = dtaset[, 5]
            x2 = dtaset[, 6]
            x3 = dtaset[, 7]
            # naive estimation
            naive_lm = lm(y ~ x1 + x2 + x3)
            naive[i, ] = naive_lm$coefficients[2:4]
            # univariate regression
            univ_lm1 = lm(y ~ x1)
            unive1[i] = univ_lm1$coefficients[2]
            univ_lm2 = lm(y ~ x2)
            unive2[i] = univ_lm2$coefficients[2]
            univ_lm3 = lm(y ~ x3)
            unive3[i] = univ_lm3$coefficients[2]
            # derived estimation
            x2dx1 = x2/x1
            x3dx1 = x3/x1
            derive_lm1 = lm(y ~ x1 + x2dx1 + x3dx1)
            derived1[i] = derive_lm1$coefficients[2]
            x1dx2 = x1/x2
            x3dx2 = x3/x2
            derive_lm2 = lm(y ~ x2 + x1dx2 + x3dx2)
            derived2[i] = derive_lm2$coefficients[2]
            x1dx3 = x1/x3
            x2dx3 = x2/x3
            derive_lm3 = lm(y ~ x3 + x1dx3 + x2dx3)
            derived3[i] = derive_lm3$coefficients[2]
          } else {
            # multivariate Cox regression
            dtaset = gend_mult_cox(n, c, trube, mut, covt, cove, beta_q)
            t = dtaset[, 1]
            status = dtaset[, 2]
            t1 = dtaset[, 3]
            t2 = dtaset[, 4]
            t3 = dtaset[, 5]
            x1 = dtaset[, 6]
            x2 = dtaset[, 7]
            x3 = dtaset[, 8]
            # naive estimation
            naive_cox = coxph(Surv(t, status) ~ x1 + x2 + x3)
            naive[i, ] = naive_cox$coefficients
            # univariate regression
            univ_cox1 = coxph(Surv(t, status) ~ x1)
            unive1[i] = univ_cox1$coefficients
            univ_cox2 = coxph(Surv(t, status) ~ x2)
            unive2[i] = univ_cox2$coefficients
            univ_cox3 = coxph(Surv(t, status) ~ x3)
            unive3[i] = univ_cox3$coefficients
            # derived estimation
            x2dx1 = x2/x1
            x3dx1 = x3/x1
            derive_cox1 = coxph(Surv(t, status) ~ x1 + x2dx1 + x3dx1)
            derived1[i] = derive_cox1$coefficients[1]
            x1dx2 = x1/x2
            x3dx2 = x3/x2
            derive_cox2 = coxph(Surv(t, status) ~ x2 + x1dx2 + x3dx2)
            derived2[i] = derive_cox2$coefficients[1]
            x1dx3 = x1/x3
            x2dx3 = x2/x3
            derive_cox3 = coxph(Surv(t, status) ~ x3 + x1dx3 + x2dx3)
            derived3[i] = derive_cox3$coefficients[1]
          }
        }
        # naive
        naive_res = colMeans(naive)
        b1e_n[j, k] = naive_res[1]
        b2e_n[j, k] = naive_res[2]
        b3e_n[j, k] = naive_res[3]
        # univariate
        b1e_un[j, k] = mean(unive1)
        b2e_un[j, k] = mean(unive2)
        b3e_un[j, k] = mean(unive3)
        # derived
        b1e_d[j, k] = mean(derived1)
        b2e_d[j, k] = mean(derived2)
        b3e_d[j, k] = mean(derived3)
        # print(c(j,k))
      }
    }
    ## summary: naive
    bias_b1_naive = b1e_n - trube[1]
    bias_b2_naive = b2e_n - trube[2]
    bias_b3_naive = b3e_n - trube[3]
    # univariate
    bias_b1_univ = b1e_un - trube[1]
    bias_b2_univ = b2e_un - trube[2]
    bias_b3_univ = b3e_un - trube[3]
    # derived
    bias_b1_deriv = b1e_d - trube[1]
    bias_b2_deriv = b2e_d - trube[2]
    bias_b3_deriv = b3e_d - trube[3]
    ####### multivariate
    bias_naive1 = cbind(rho_e_set, bias_b1_naive)
    bias_naive2 = cbind(rho_e_set, bias_b2_naive)
    bias_naive3 = cbind(rho_e_set, bias_b3_naive)
    colnames(bias_naive1) = colnames(bias_naive2) = colnames(bias_naive3) = c("rho_e",
                                                                              paste0("rho_T=", seq(0.1, 0.4, by = 0.1)))
    rownames(bias_naive1) = rownames(bias_naive2) = rownames(bias_naive3) = rep("Multivariate",
                                                                                4)
    b1_biva = rbind(b1_biva, bias_naive1)
    b2_biva = rbind(b2_biva, bias_naive2)
    b3_biva = rbind(b3_biva, bias_naive3)
    # univariate
    bias_univ1 = cbind(rho_e_set, bias_b1_univ)
    bias_univ2 = cbind(rho_e_set, bias_b2_univ)
    bias_univ3 = cbind(rho_e_set, bias_b3_univ)
    colnames(bias_univ1) = colnames(bias_univ2) = colnames(bias_univ3) = c("rho_e", paste0("rho_T=",
                                                                                           seq(0.1, 0.4, by = 0.1)))
    rownames(bias_univ1) = rownames(bias_univ2) = rownames(bias_univ3) = rep("Univariate",
                                                                             4)
    b1_uni = rbind(b1_uni, bias_univ1)
    b2_uni = rbind(b2_uni, bias_univ2)
    b3_uni = rbind(b3_uni, bias_univ3)
    # derived
    bias_deriv1 = cbind(rho_e_set, bias_b1_deriv)
    bias_deriv2 = cbind(rho_e_set, bias_b2_deriv)
    bias_deriv3 = cbind(rho_e_set, bias_b3_deriv)
    colnames(bias_deriv1) = colnames(bias_deriv2) = colnames(bias_deriv3) = c("rho_e",
                                                                              paste0("rho_T=", seq(0.1, 0.4, by = 0.1)))
    rownames(bias_deriv1) = rownames(bias_deriv2) = rownames(bias_deriv3) = rep("Derived",
                                                                                4)
    b1_derive = rbind(b1_derive, bias_deriv1)
    b2_derive = rbind(b2_derive, bias_deriv2)
    b3_derive = rbind(b3_derive, bias_deriv3)
  }
  
  # univariate estimation results for each \rho_{\epsilon} are almost the same
  b1_final_res = rbind(b1_biva[1:4, ], b1_uni[2, ], b1_derive[1:4, ], 
  b1_biva[5:8, ], b1_uni[6, ], b1_derive[5:8, ])
  b2_final_res = rbind(b2_biva[1:4, ], b2_uni[2, ], b2_derive[1:4, ], 
  b2_biva[5:8, ], b2_uni[6, ], b2_derive[5:8, ])
  b3_final_res = rbind(b3_biva[1:4, ], b3_uni[2, ], b3_derive[1:4, ], 
  b3_biva[5:8, ], b3_uni[6, ], b3_derive[5:8, ])
  # estimation bias of \beta_1
  res = cbind(b1_final_res, b2_final_res[, -1], b3_final_res[, -1])
  return(res)
}

# bivariate linear regression for observations with zero values, this program can
# reproduce results in Tables S6-S7 of Web Appendix 
# parameter description 
# trube: true regression coefficients 
# sigma_x1: standard error of observed covariate x1 
# sigma_x2: standard error of observed covariate x2 
# mis_p: proportion of zero values
est_biva_linear_with_zero_obs <- function(trube, sigma_x1, sigma_x2, mis_p) {
  scenario_set = 1:4
  m = n * mis_p  #number of observations with zero values
  b1_biva = b2_biva = b1_uni = b2_uni = b1_derive = b2_derive = b1_derivef = b2_derivef = NULL
  for (scenario in scenario_set) {
    if (scenario == 1) {
      mut = c(1, 1)
    } else if (scenario == 2) {
      mut = c(5, 5)
    } else if (scenario == 3) {
      mut = c(5, 10)
    } else if (scenario == 4) {
      mut = c(10, 5)
    }
    
    sigma_t1 = sigma_x1 * rho_1  #square-root of variance of true variable X1
    sigma_t2 = sigma_x2 * rho_2  #square-root of variance of true variable X2
    sigma_e1 = sigma_x1 * sqrt(1 - rho_1^2)  #square-root of variance of error variable e1
    sigma_e2 = sigma_x2 * sqrt(1 - rho_2^2)  #square-root of variance of error variable e2
    
    # set of \rho_T
    rho_t_set = c(0.3, 0.5, 0.7, 0.9)
    len_rhot = length(rho_t_set)
    # set of \rho_{\epsilon}
    rho_e_set = seq(0.5, 0.9, by = 0.1)
    len_rhoe = length(rho_e_set)
    b1e_n = b2e_n = b1e_un = b2e_un = b1e_d = b2e_d = b1e_df = b2e_df = matrix(0, len_rhoe,
                                                                               len_rhot)
    for (k in 1:len_rhot) {
      rho_t = rho_t_set[k]
      for (j in 1:len_rhoe) {
        rho_e = rho_e_set[j]
        covt = matrix(c(sigma_t1^2, rho_t * sigma_t1 * sigma_t2, rho_t * sigma_t1 *
                          sigma_t2, sigma_t2^2), 2, 2)
        cove = matrix(c(sigma_e1^2, rho_e * sigma_e1 * sigma_e2, rho_e * sigma_e1 *
                          sigma_e2, sigma_e2^2), 2, 2)
        naive = matrix(0, N, ncol = 2)
        unive1 = unive2 = derived1 = derived2 = derived1f = derived2f = rep(N, 0)
        for (i in 1:N) {
          set.seed(123456 + i)
          dtaset = gend_biva_linear(n, trube, mut, covt, cove)
          y = dtaset[, 1]
          t1 = dtaset[, 2]
          t2 = dtaset[, 3]  #true covariates
          x1 = dtaset[, 4]
          x2 = dtaset[, 5]  #observed covariates
          index = 1:m  # the first m observations can be treated a random subset of the n sample size
          x1[index] = 0
          x2[index] = 0
          # naive estimation
          naive_lm = lm(y ~ x1 + x2)
          naive[i, ] = naive_lm$coefficients[2:3]
          # univariate regression
          univ_lm1 = lm(y ~ x1)
          unive1[i] = univ_lm1$coefficients[2]
          univ_lm2 = lm(y ~ x2)
          unive2[i] = univ_lm2$coefficients[2]
          # derived1 estimation (first strategy, i.e., replacing undefined
          # derived variable Z by zero
          z1 = numeric(n)
          z1[(m + 1):n] = x2[(m + 1):n]/x1[(m + 1):n]
          derive_lm1f = lm(y ~ x1 + z1)
          derived1f[i] = derive_lm1f$coefficients[2]
          z2 = numeric(n)
          z2[(m + 1):n] = x1[(m + 1):n]/x2[(m + 1):n]
          derive_lm2f = lm(y ~ z2 + x2)
          derived2f[i] = derive_lm2f$coefficients[3]
          # derived2 (the second strategy, i.e., deleting observations with
          # zero values)
          dtaset0 = dtaset[-index, ]
          x1d = dtaset0[, 4]
          x2d = dtaset0[, 5]
          yd = dtaset0[, 1]
          x2dx1 = x2d/x1d
          derive_lm1 = lm(yd ~ x1d + x2dx1)
          derived1[i] = derive_lm1$coefficients[2]
          x1dx2 = x1d/x2d
          derive_lm2 = lm(yd ~ x1dx2 + x2d)
          derived2[i] = derive_lm2$coefficients[3]
          
        }
        naive_res = colMeans(naive)
        b1e_n[j, k] = naive_res[1]
        b2e_n[j, k] = naive_res[2]
        b1e_un[j, k] = mean(unive1)
        b2e_un[j, k] = mean(unive2)
        b1e_df[j, k] = mean(derived1f)
        b2e_df[j, k] = mean(derived2f)
        b1e_d[j, k] = mean(derived1)
        b2e_d[j, k] = mean(derived2)
        print(c(j, k))
      }
    }
    ## summary simulation results
    bias_b1_naive = b1e_n - trube[1]
    bias_b2_naive = b2e_n - trube[2]
    bias_b1_univ = b1e_un - trube[1]
    bias_b2_univ = b2e_un - trube[2]
    bias_b1_derivf = b1e_df - trube[1]
    bias_b2_derivf = b2e_df - trube[2]
    bias_b1_deriv = b1e_d - trube[1]
    bias_b2_deriv = b2e_d - trube[2]
    
    # naive
    bias_naive1 = cbind(rho_e_set, bias_b1_naive)
    bias_naive2 = cbind(rho_e_set, bias_b2_naive)
    colnames(bias_naive1) = colnames(bias_naive2) = c("rho_e", paste0("rho_T=", seq(0.3,
                                                                                    0.9, by = 0.2)))
    rownames(bias_naive1) = rownames(bias_naive2) = rep("Bivariate", 5)
    b1_biva = rbind(b1_biva, bias_naive1)
    b2_biva = rbind(b2_biva, bias_naive2)
    # univariate
    bias_univ1 = cbind(rho_e_set, bias_b1_univ)
    bias_univ2 = cbind(rho_e_set, bias_b2_univ)
    colnames(bias_univ1) = colnames(bias_univ2) = c("rho_e", paste0("rho_T=", seq(0.3,
                                                                                  0.9, by = 0.2)))
    rownames(bias_univ1) = rownames(bias_univ2) = rep("Univariate", 5)
    b1_uni = rbind(b1_uni, bias_univ1)
    b2_uni = rbind(b2_uni, bias_univ2)
    # derived1
    bias_deriv1f = cbind(rho_e_set, bias_b1_derivf)
    bias_deriv2f = cbind(rho_e_set, bias_b2_derivf)
    colnames(bias_deriv1f) = colnames(bias_deriv2f) = c("rho_e", paste0("rho_T=", seq(0.3,
                                                                                      0.9, by = 0.2)))
    rownames(bias_deriv1f) = rownames(bias_deriv2f) = rep("Derived1", 5)
    b1_derivef = rbind(b1_derivef, bias_deriv1f)
    b2_derivef = rbind(b2_derivef, bias_deriv2f)
    # derived2
    bias_deriv1 = cbind(rho_e_set, bias_b1_deriv)
    bias_deriv2 = cbind(rho_e_set, bias_b2_deriv)
    colnames(bias_deriv1) = colnames(bias_deriv2) = c("rho_e", paste0("rho_T=", seq(0.3,
                                                                                    0.9, by = 0.2)))
    rownames(bias_deriv1) = rownames(bias_deriv2) = rep("Derived2", 5)
    b1_derive = rbind(b1_derive, bias_deriv1)
    b2_derive = rbind(b2_derive, bias_deriv2)
  }
  
  # Results in Tables S6
  b1_res_tables6 = rbind(b1_biva[1:5, ], Univariate = b1_uni[3, ], b1_derivef[1:5, ], 
  b1_derive[1:5, ], b1_biva[6:10, ], Univariate = b1_uni[8, ], b1_derivef[6:10, ],
  b1_derive[6:10, ])
  b2_res_tables6 = rbind(b2_biva[1:5, ], Univariate = b2_uni[3, ], b2_derivef[1:5, ], 
  b2_derive[1:5, ], b2_biva[6:10, ], Univariate = b2_uni[8, ], b2_derivef[6:10, ], 
  b2_derive[6:10, ])
  tabs6 = cbind(b1_res_tables6, b2_res_tables6[, -1])
  
  # Results in Table S7
  b1_res_tables7 = rbind(b1_biva[11:15, ], Univariate = b1_uni[13, ], 
  b1_derivef[11:15, ], b1_derive[11:15, ], b1_biva[16:20, ], 
  Univariate = b1_uni[18, ], b1_derivef[16:20, ], b1_derive[16:20, ])
  b2_res_tables7 = rbind(b2_biva[11:15, ], Univariate = b2_uni[13, ], 
  b2_derivef[11:15, ], b2_derive[11:15, ], b2_biva[16:20, ], 
  Univariate = b2_uni[18, ], b2_derivef[16:20, ], b2_derive[16:20, ])
  tabs7 = cbind(b1_res_tables7, b2_res_tables7[, -1])
  
  res = list(tables6 = tabs6, tables7 = tabs7)
  return(res)
}

# Calculate Concordance rate, bias and MSE of four methods of estimation, which can
# reproduce results in Table 7 of manuscript
cal_cr_bias_mse <- function(scen_set) {
  final_res = NULL
  for (scenario in scen_set) {
    # scenario = 1 #scenario=1 is for male; scenario=2 is for female X1 = protein,
    # x2 = fat, X3 = energy, male
    if (scenario == 1) {
      mut = c(9.356683, 8.654661, 11.937024)
      trube = c(0.172, 0.228, 0.187)  #
      covt = matrix(c(1.05927134, 0.050755983, 0.33722698, 0.05075598, 0.676466609,
                      0.31167934, 0.33722698, 0.31167934, 1.1386985), byrow = T, 3, 3)
      covt = round(covt, 4)
      cove = matrix(c(0.8575957, 0.7162641, 0.7517087, 0.7162641, 0.8836357, 0.7924088,
                      0.7517087, 0.7924088, 0.9131195), byrow = T, 3, 3)
      cove = round(cove, 4)
      sigtx = matrix(c(0.38838743, 0.021051055, 0.093149182, 0.01860995, 0.280564669,
                       0.086092387, 0.12364605, 0.129269071, 0.31453247), byrow = T, 3, 3)
      sigtx = round(sigtx, 4)
      beta_q = c(0.3667, 0.4148, 0.2762)
    } else if (scenario == 2) {
      # female
      mut = c(9.245602, 7.510393, 11.167943)
      trube = c(0.16, 0.112, 0.074)
      covt = matrix(c(1.07312719, 0.05953603, 0.32754194, 0.05953603, 1.07539124, 0.4410434,
                      0.32754194, 0.4410434, 1.17795383), byrow = T, 3, 3)
      covt = round(covt, 4)
      cove = matrix(c(0.9046743, 0.7202257, 0.7603176, 0.7202257, 0.9308771, 0.8315039,
                      0.7603176, 0.8315039, 0.9254783), byrow = T, 3, 3)
      cove = round(cove, 4)
      sigtx = matrix(c(0.31983837, 0.015094115, 0.082384259, 0.01774431, 0.272642958,
                       0.110932463, 0.09762168, 0.11181733, 0.296282223), byrow = T, 3, 3)
      sigtx = round(sigtx, 4)
      beta_q = c(0.298, 0.2535, 0.2515)
    }
    
    if (scenario < 2) {
      n = 5785  #sample size for male
    } else {
      n = 9472  #sample size for female
    }
    ## keep simulation results
    naive = naive_p = univ = univ_p = derived = derived_p = rce = rce_p = matrix(0, N,
                                                                                 ncol = 3)
    
    # iteration
    for (i in 1:N) {
      set.seed(10000 + i)
      dtaset = gend_mult_cox(n, c = 1.6, trube, mut, covt, cove, beta_q)
      t = dtaset[, 1]
      delta = dtaset[, 2]
      t1 = dtaset[, 3]
      t2 = dtaset[, 4]
      t3 = dtaset[, 5]
      x1 = dtaset[, 6]
      x2 = dtaset[, 7]
      x3 = dtaset[, 8]
      # naive estimation
      naive_cox = coxph(Surv(t, delta) ~ x1 + x2 + x3)
      naive_summ = summary(naive_cox)$coefficients
      naive[i, ] = naive_summ[, 1]
      naive_p[i, ] = naive_summ[, 5]
      # rc estimation
      x = cbind(x1, x2, x3)
      sigxx = cov(x)  
      umx = apply(x, 2, mean)
      etx = matrix(rep(mut, n), byrow = TRUE, ncol = 3) + t(sigtx %*% solve(sigxx) %*%
                                                              t(x - umx))
      m1 = etx[, 1]
      m2 = etx[, 2]
      m3 = etx[, 3]
      rc_cox = coxph(Surv(t, delta) ~ m1 + m2 + m3)
      rc_summ = summary(rc_cox)$coefficients
      rce[i, ] = rc_summ[, 1]
      rce_p[i, ] = rc_summ[, 5]
      # univariate estimation
      univ_cox1 = coxph(Surv(t, delta) ~ x1)
      univ_summ1 = summary(univ_cox1)$coefficients
      univ[i, 1] = univ_summ1[1, 1]
      univ_p[i, 1] = univ_summ1[1, 5]
      
      univ_cox2 = coxph(Surv(t, delta) ~ x2)
      univ_summ2 = summary(univ_cox2)$coefficients
      univ[i, 2] = univ_summ2[1, 1]
      univ_p[i, 2] = univ_summ2[1, 5]
      
      univ_cox3 = coxph(Surv(t, delta) ~ x3)
      univ_summ3 = summary(univ_cox3)$coefficients
      univ[i, 3] = univ_summ3[1, 1]
      univ_p[i, 3] = univ_summ3[1, 5]
      # derived estimation
      x1d = x1/x3
      x2d = x2/x3
      derive_cox3 = coxph(Surv(t, delta) ~ x1d + x2d + x3)
      derive_summ3 = summary(derive_cox3)$coefficients
      derived[i, 3] = derive_summ3[3, 1]
      derived_p[i, 3] = derive_summ3[3, 5]
      
      x1df = x1/x2
      x3df = x3/x2
      derive_cox2 = coxph(Surv(t, delta) ~ x1df + x2 + x3df)
      derive_summ2 = summary(derive_cox2)$coefficients
      derived[i, 2] = derive_summ2[2, 1]
      derived_p[i, 2] = derive_summ2[2, 5]
      
      x2dff = x2/x1
      x3dff = x3/x1
      derive_cox1 = coxph(Surv(t, delta) ~ x1 + x2dff + x3dff)
      derive_summ1 = summary(derive_cox1)$coefficients
      derived[i, 1] = derive_summ1[1, 1]
      derived_p[i, 1] = derive_summ1[1, 5]
    }
    
    ### summary estimation bias
    # naive
    naive_res = colMeans(naive)
    naive_bias = naive_res - trube
    # univariate
    univ_res = colMeans(univ)
    univ_bias = univ_res - trube
    # derived
    derive_res = colMeans(derived)
    derived_bias = derive_res - trube
    # rc
    rc_res = colMeans(rce)
    rc_bias = rc_res - trube
    Bias_res = data.frame(Naive = naive_bias, Univariate = univ_bias, Derived = derived_bias,
                          RC = rc_bias)
    row.names(Bias_res) = paste0("beta", 1:3)
    
    ##### MSE
    trubem = matrix(rep(trube, N), byrow = T, ncol = 3)
    naive_mse = colMeans((naive - trubem)^2)
    univ_mse = colMeans((univ - trubem)^2)
    derived_mse = colMeans((derived - trubem)^2)
    rc_mse = colMeans((rce - trubem)^2)
    MSE_res = data.frame(Naive = naive_mse, Univariate = univ_mse, Derived = derived_mse,
                         RC = rc_mse)
    row.names(MSE_res) = paste0("beta", 1:3)
    
    ##### Concordance rate
    cal_cr = function(est_p, rce_p) {
      # est_p can be naive_p, univ_p and derived_p, which is a N*3 matrix
      cr_res = numeric(3)
      for (j in 1:3) {
        # comparing results with RC
        n00 = sum(est_p[, j] >= 0.05 & rce_p[, j] >= 0.05)  #both cannot reject
        n01 = sum(est_p[, j] >= 0.05 & rce_p[, j] < 0.05)  #naive cannot reject but rc can reject
        n10 = sum(est_p[, j] < 0.05 & rce_p[, j] >= 0.05)  #naive reject but rc cannot reject
        n11 = sum(est_p[, j] < 0.05 & rce_p[, j] < 0.05)  #both reject
        cr_res[j] = (n00 + n11)/N  #ratio of both reject and both cannot reject
      }
      return(cr_res)
    }
    
    naive_cr = cal_cr(naive_p, rce_p)
    univ_cr = cal_cr(univ_p, rce_p)
    derived_cr = cal_cr(derived_p, rce_p)
    
    CR_res = data.frame(Naive = naive_cr, Univariate = univ_cr, Derived = derived_cr)
    row.names(CR_res) = paste0("beta", 1:3)
    
    # keep results
    temp_res = cbind(CR_res, Bias_res, MSE_res)
    final_res = rbind(final_res, temp_res)
  }
  colnames(final_res) = c("multivariate", "Univariate", "Derived", rep(c("multivariate",
                                                                         "Univariate", "Derived", "RC"), 2))
  return(final_res)
}
