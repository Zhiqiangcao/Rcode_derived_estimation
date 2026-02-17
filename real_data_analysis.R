# Real data analysis 

library(survival)
library(rootSolve)

# input source code for real data analysis
source("./real_data_analysis_functions.R")

# Table 5 load data
load("./epic_mask.RData")

est_tab5 = four_methods_est(epic_mask)
saveRDS(est_tab5, "./results/est_tab5.rds")
est_tab5 = readRDS("./results/est_tab5.rds")
colnames(est_tab5) = c(paste0("Multivariate_", c("Est", "SE", "P_value")), paste0("Univariate_",
                                                                                  c("Est", "SE", "P_value")), paste0("Derived_", c("Est", "SE", "P_value")), paste0("RC_",
                                                                                                                                                                    c("Est", "SE", "P_value")))
row.names(est_tab5) = c(paste0("Male_", c("Protein", "Fat", "Energy", "BMI", "Edu", "Smoke")),
                        paste0("Female_", c("Protein", "Fat", "Energy", "BMI", "Edu", "Smoke")))
write.csv(est_tab5, file = "./results/Table5.csv")


# Table 6 load data
load("./validation_mask.RData")
# male's validation study result
male_res = valid_resu(validation_mask, sex = 1, rho = 0.7, bhi = 0.6, bhj = 0.6)
saveRDS(male_res, "./results/male_res.rds")
male_res = readRDS("./results/male_res.rds")
male_b = round(male_res$b, 3)  #vector b in Table 6 for male
male_mut = round(male_res$mut, 3)  #vector mu_t in Table 6 for male
male_sigma_t = male_res$Sigma_T[2:4, 2:4]  #estimated variance-covariance of true variables
# compute correlation-coefficients in off-diagonal elements
male_sigma_T = cal_corr_matrix(male_sigma_t)
male_sigma_T = round(male_sigma_T, 3)
male_sigma_e = male_res$Sigma_e[2:4, 2:4]  #estimated variance-covariance of error variables
male_sigma_eps = cal_corr_matrix(male_sigma_e)
male_sigma_eps = round(male_sigma_eps, 3)


# female's validation study result
female_res = valid_resu(validation_mask, sex = 2, rho = 0.7, bhi = 0.6, bhj = 0.6)
saveRDS(female_res, "./results/female_res.rds")
female_res = readRDS("./results/female_res.rds")
female_b = round(female_res$b, 3)  #vector b in Table 6 for female
female_mut = round(female_res$mut, 3)  #vector mu_t in Table 6 for female
female_sigma_t = female_res$Sigma_T[2:4, 2:4]  #estimated variance-covariance of true variables
# compute correlation-coefficients in off-diagonal elements
female_sigma_T = cal_corr_matrix(female_sigma_t)
female_sigma_T = round(female_sigma_T, 3)
female_sigma_e = female_res$Sigma_e[2:4, 2:4]  #estimated variance-covariance of error variables
female_sigma_eps = cal_corr_matrix(female_sigma_e)
female_sigma_eps = round(female_sigma_eps, 3)

# note: \beta in Table 6 are RC estimates of protein, fat and energy in Table 5
male_beta = est_tab5[1:3, 10]
female_beta = est_tab5[7:9, 10]

Table6 = list(male_b = male_b, female_b = female_b, male_beta = male_beta, female_beta = female_beta,
              male_mut = male_mut, female_mut = female_mut, male_sigma_T = male_sigma_T, female_sigma_T = female_sigma_T,
              male_sigma_eps = male_sigma_eps, female_sigma_eps = female_sigma_eps)
write.csv(Table6, file = "./results/Table6.csv")

# Note: since we use masked main and validation study data, results of Tables 5-6 are 
# slgihtly different from those in paper 