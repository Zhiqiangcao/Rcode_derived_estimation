### simulation studies:

library(mvtnorm)
library(survival)
library(glmnet)


# input source code for simulation studies
source("./simulation_functions.R")

# parameter settings for all simulations
N = 1000  #simulation times
n = 400  #sample size for each simulation 
rho_1 = 0.5  #reliability index of X1
rho_2 = 0.6  #reliability index of X2

# row and column names used in corresponding Tables
colname1 = c("rho_e", paste0("beta1_rho_T=", seq(0.3, 0.9, by = 0.2)), paste0("beta2_rho_T=",
                                                                              seq(0.3, 0.9, by = 0.2)))
colname2 = c("rho_e", paste0("beta1_rho_T=", seq(0.1, 0.4, by = 0.1)), paste0("beta2_rho_T=",
                                                                              seq(0.1, 0.4, by = 0.1)), paste0("beta3_rho_T=", seq(0.1, 0.4, by = 0.1)))
rowname1 = c(rep("Bivariate", 5), rep("Ridge", 5), "Univariate", rep("mu1=mu2=1", 5), rep("mu1=mu2=5",
                                                                                          5), rep("mu1=5_mu2=10", 5), rep("mu1=10_mu2=5", 5), rep("mu1=mu2=10", 5))
rowname2 = c(rep("Bivariate", 5), "Univariate", rep("mu1=mu2=1", 5), rep("mu1=mu2=5", 5),
             rep("mu1=5_mu2=10", 5), rep("mu1=10_mu2=5", 5), rep("mu1=mu2=10", 5))
rowname3 = c(rep("Male_Multivariate", 4), "Male_Univariate", rep("Male_Derived", 4), rep("Female_Multivariate",
                                                                                         4), "Female_Univariate", rep("Female_Derived", 4))

# Table 2
est_tab2 = est_biva_linear(trube = c(1, 2), sigma_x1 = 1, sigma_x2 = 1, error_type = 1, ridge_reg = TRUE)
saveRDS(est_tab2, "./results/est_tab2.rds")
# note it will take about 1.26 hour
est_tab2 = readRDS("./results/est_tab2.rds")
colnames(est_tab2) = colname1
row.names(est_tab2) = rowname1
write.csv(round(est_tab2, 3), file = "./results/Table2.csv")


# Table 3
est_tab3 = est_biva_linear(trube = c(2, 1), sigma_x1 = 1, sigma_x2 = 1, error_type = 1, ridge_reg = TRUE)
saveRDS(est_tab3, "./results/est_tab3.rds")
# note it will take about 1.36 hour
est_tab3 = readRDS("./results/est_tab3.rds")
colnames(est_tab3) = colname1
row.names(est_tab3) = rowname1
write.csv(round(est_tab3, 3), file = "./results/Table3.csv")


# Table 4
est_tab4 = est_biva_cox(trube = c(1, 2), sigma_x1 = 1, sigma_x2 = 1)
saveRDS(est_tab4, "./results/est_tab4.rds")
# note it will take about 5.9 hour
est_tab4 = readRDS("./results/est_tab4.rds")
colnames(est_tab4) = colname1
row.names(est_tab4) = c(rep("mu1=mu2=1_Bivariate", 5), rep("mu1=mu2=1_Ridge", 5), "mu1=mu2=1_Univariate",
                        rep("mu1=mu2=1_Derived", 5), rep("mu1=mu2=5_Bivariate", 5), rep("mu1=mu2=5_Ridge", 5),
                        "mu1=mu2=5_Univariate", rep("mu1=mu2=5_Derived", 5))
write.csv(round(est_tab4, 3), file = "./results/Table4.csv")

# Table 7
est_tab7 = cal_cr_bias_mse(scen_set = 1:2)
saveRDS(est_tab7, "./results/est_tab7.rds")
# note it will take about 11 minutes
est_tab7 = readRDS("./results/est_tab7.rds")
colnames(est_tab7) = c(paste0("CR_", c("Multivariate", "Univariate", "Derived")), paste0("Bias_",
                                                                                         c("Multivariate", "Univariate", "Derived", "RC")), paste0("MSE_", c("Multivariate", "Univariate",
                                                                                                                                                             "Derived", "RC")))
row.names(est_tab7) = c(paste0("Male_", paste0("beta", 1:3)), paste0("Female_", paste0("beta",
                                                                                       1:3)))
write.csv(round(est_tab7, 3), file = "./results/Table7.csv")


# Table S1
est_tabs1 = est_biva_linear(trube = c(1, -2), sigma_x1 = 1, sigma_x2 = 1, error_type = 1,
                            ridge_reg = FALSE)
saveRDS(est_tabs1, "./results/est_tabs1.rds")
# note it will take about 7 minutes
est_tabs1 = readRDS("./results/est_tabs1.rds")
colnames(est_tabs1) = colname1
row.names(est_tabs1) = rowname2
write.csv(round(est_tabs1, 3), file = "./results/Tables1.csv")

# Table S2
est_tabs2 = est_biva_linear(trube = c(1, 1.5), sigma_x1 = 1, sigma_x2 = 1, error_type = 1,
                            ridge_reg = FALSE)
saveRDS(est_tabs2, "./results/est_tabs2.rds")
# note it will take about 7 minutes
est_tabs2 = readRDS("./results/est_tabs2.rds")
colnames(est_tabs2) = colname1
row.names(est_tabs2) = rowname2
write.csv(round(est_tabs2, 3), file = "./results/Tables2.csv")


# Table S3
est_tabs3 = est_biva_linear(trube = c(1, 2), sigma_x1 = sqrt(10), sigma_x2 = sqrt(10), error_type = 1,
                            ridge_reg = FALSE)
saveRDS(est_tabs3, "./results/est_tabs3.rds")
# note it will take about 7 minutes
est_tabs3 = readRDS("./results/est_tabs3.rds")
colnames(est_tabs3) = colname1
row.names(est_tabs3) = rowname2
write.csv(round(est_tabs3, 3), file = "./results/Tables3.csv")


# Table S4
est_tabs4 = est_mult_regression(linear = TRUE)
saveRDS(est_tabs4, "./results/est_tabs4.rds")
# note it will take about 3 minutes
est_tabs4 = readRDS("./results/est_tabs4.rds")
colnames(est_tabs4) = colname2
row.names(est_tabs4) = rowname3
write.csv(round(est_tabs4, 3), file = "./results/Tables4.csv")


# Table S5
est_tabs5 = est_mult_regression(linear = FALSE)
saveRDS(est_tabs5, "./results/est_tabs5.rds")
# note it will take about 12 minutes
est_tabs5 = readRDS("./results/est_tabs5.rds")
colnames(est_tabs5) = colname2
row.names(est_tabs5) = rowname3
write.csv(round(est_tabs5, 3), file = "./results/Tables5.csv")

# Table S6-S7
est_tabs6_s7 = est_biva_linear_with_zero_obs(trube = c(1, 2), sigma_x1 = 1, sigma_x2 = 1,
                                             mis_p = 0.1)
saveRDS(est_tabs6_s7, "./results/est_tabs6_s7.rds")
# note it will take about 6 minutes
est_tabs6_s7 = readRDS(".results/est_tabs6_s7.rds")
est_tabs6 = est_tabs6_s7$tables6
colnames(est_tabs6) = colname1
row.names(est_tabs6) = c(rep("mu1=mu2=1_Bivariate", 5), "mu1=mu2=1_Univariate", rep("mu1=mu2=1_Derived1",
                                                                                    5), rep("mu1=mu2=1_Derived2", 5), rep("mu1=mu2=5_Bivariate", 5), "mu1=mu2=5_Univariate",
                         rep("mu1=mu2=5_Derived1", 5), rep("mu1=mu2=5_Derived2", 5))
write.csv(round(est_tabs6, 3), file = ".results/Tables6.csv")

est_tabs7 = est_tabs6_s7$tables7
colnames(est_tabs7) = colname1
row.names(est_tabs7) = c(rep("mu1=5_mu2=10_Bivariate", 5), "mu1=5_mu2=10_Univariate", rep("mu1=5_mu2=10_Derived1",
                                                                                          5), rep("mu1=5_mu2=10_Derived2", 5), rep("mu1=10_mu2=5_Bivariate", 5), "mu1=10_mu2=5_Univariate",
                         rep("mu1=10_mu2=5_Derived1", 5), rep("mu1=10_mu2=5_Derived2", 5))
write.csv(round(est_tabs7, 3), file = "./results/Tables7.csv")

# Table S8
est_tabs8 = est_biva_linear(trube = c(1, 2), sigma_x1 = 1, sigma_x2 = 1, error_type = 2,
                            ridge_reg = FALSE)
saveRDS(est_tabs8, "./results/est_tabs8.rds")
# note it will take about 2 minutes
est_tabs8 = readRDS(".results/est_tabs8.rds")
colnames(est_tabs8) = colname1
row.names(est_tabs8) = c(rep("mu1=mu2=1_Bivariate", 5), "mu1=mu2=1_Univariate", rep("mu1=mu2=1_Derived",
                                                                                    5), rep("mu1=mu2=5_Bivariate", 5), "mu1=mu2=5_Univariate", rep("mu1=mu2=5_Derived", 5))
write.csv(round(est_tabs8, 3), file = ".results/Tables8.csv")