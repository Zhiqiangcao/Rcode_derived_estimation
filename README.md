# Rcode_derived_estimation
The Rcode is to reproduce simulation results in "Reducing the effect of correlated additive measurement errors on regression models with no ancillary study" By Zhiqiang Cao and Man Yu Wong.


For questions or comments about the code, please contact Zhiqiang Cao zcaoae@connect.ust.hk.  
This folder includes the following functions:

The numerical study has two parts: 

Part (1) simulation

R scripts: main_simulations.R, simulation_functions.R  

Run script main_simulations.R to produce Tables 2-4, 7 in manuscript and Tables S1-S8 in Web Appendix

Script simulation_functions.R defines the functions used in the simulation studies.

(The result saved as "est_tabi.rds" (i=2,3,4,7) or "est_tabsj.rds" (j=1,2,3,4,5,6,7,8), which can be used to generate  Tables 2-4,7 and Tables S1-S8)


Part (2) data example

R scripts: real_data_analysis.R, real_data_analysis_functions.R and epic_mask.R

Data set:  epic_mask.Rdata and validation_mask.Rdata (The details of the mask procedure are provided in R script epic_mask.R)

Run real_data_analysis.R to produce Tables 5-6 in manuscript

Script real_data_analysis_functions.R defines the functions used in the real data analysis.

(The result saved as "est_tab5.rds", "male_res.rds" and "female_res.rds", which can be used to generate Tables 5-6)



Note (1): Tables and Figures are saved in the results folder.

Note (2):  You will see warning messages when producing Table 4 and Table 7, the warning message may like this

"In coxph.fit(X, Y, istrat, offset, init, control, weights = weights,  ... : Loglik converged before variable  2 ; coefficient may be infinite." 

Note (3): Figure 1 of this paper is quoted from Figure 1(a) in the reference "Day, N.E., Wong, M.Y., Bingham, S., Khaw, K.T., Luben, 

R., Michels, K.B., Welch, A., & Wareham, N.J.  (2004). Correlated measurement error-implications for nutritional epidemiology. Int J Epidemiol, 33, 1373-1381."




