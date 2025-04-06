#Real data analysis based on multivariate, univariate, RC and derived estimation in paper.
###note:  real data has been masked (they are different to real data used in paper) since data is confidential
rm(list=ls()) #clean environment
library(MASS)
library(mvtnorm)
library(survival)
library(nleqslv)

#make sure the current working directory is the code folder, i.e., 
#change directory where unzipped files are located

load("epic_mask.RData")

setwd("C:/Users/82655/Dropbox/research/measurment_error/paper1_amle_logistic/paper1_supplementary_Code_and_Data")

#analyze male and female separately 
for (sex in 1:2){
  if(sex == 1) epic = epic_mask[epic_mask$sex==1,]  #male
  if(sex == 2) epic = epic_mask[epic_mask$sex==2,]  #female
  
  #sqrt-root transformation, then divided by their standard deviation
  prot = sqrt(epic$qe_prot)
  prot = prot/sd(prot)
  fat = sqrt(epic$qe_fat)
  fat = fat/sd(fat)
  energy = sqrt(epic$qe_energy)
  energy = energy/sd(energy)
  
  ##error-free covariates setup
  bmi = scale(epic$bmi_adj); edu = 1*(epic$l_school>=3)
  smoke = 1*(epic$smoke_stat==3) 
  
  #center
  is11=epic$centre==11; is12=epic$centre==12; is13=epic$centre==13
  is14=epic$centre==14; is15=epic$centre==15; is16=epic$centre==16
  is21=epic$centre==21; is22=epic$centre==22; is23=epic$centre==23
  is24=epic$centre==24; is25=epic$centre==25;
  is31=epic$centre==31; is32=epic$centre==32; is33=epic$centre==33
  is34=epic$centre==34; is35=epic$centre==35;
  is41=epic$centre==41; is42=epic$centre==42; is51=epic$centre==51
  is52=epic$centre==52; is71=epic$centre==71; is72=epic$centre==72
  
  ###construct E(x|w) and 
  w = cbind(prot,fat,energy)
  sigww = cov(w)  #variance-covariance between variables measured with error
  umw = apply(w,2,mean)
  n = dim(w)[1]
  
  ###########################naive estimation
  naive1 = coxph(Surv(st,event)~prot+fat+energy+bmi+edu+smoke+factor(centre),data=epic)
  est1 = summary(naive1);
  naive1_res = round(est1$coefficients[1:6,-2],3); 
  
  ###########################RC estimation
  if(sex==1){
    sigxx=matrix(c(1.05927134,0.050755983,0.33722698,
                   0.05075598,0.676466609,0.31167934,
                   0.33722698,0.31167934,1.13869850),byrow=T,3,3)
    sigxx = round(sigxx,4)
    siguu = matrix(c(0.8575957,0.7162641,0.7517087,
                     0.7162641,0.8836357,0.7924088,
                     0.7517087,0.7924088,0.9131195),byrow=T,3,3)
    siguu = round(siguu,4)
    sigxw = matrix(c(0.38838743,0.021051055,0.093149182,
                     0.01860995,0.280564669,0.086092387,
                     0.12364605,0.129269071,0.314532470),byrow=T,3,3)
    sigxw = round(sigxw,4)
    beta=c(0.3667,0.4148,0.2762)
    mux = c(9.356683,8.654661,11.937024)
  }else if(sex==2){
    sigxx=matrix(c(1.07312719,0.05953603,0.32754194,
                   0.05953603,1.07539124,0.44104340,
                   0.32754194,0.44104340,1.17795383),byrow=T,3,3)
    sigxx = round(sigxx,4)
    siguu=matrix(c(0.9046743,0.7202257,0.7603176,
                   0.7202257,0.9308771,0.8315039,
                   0.7603176,0.8315039,0.9254783),byrow=T,3,3)
    siguu = round(siguu,4)
    sigxw=matrix(c(0.31983837,0.015094115,0.082384259,
                   0.01774431,0.272642958,0.110932463,
                   0.09762168,0.111817330,0.296282223),byrow=T,3,3)
    sigxw = round(sigxw,4)
    beta=c(0.2980,0.2535,0.2515)
    mux = c(9.245602,7.510393,11.167943)
  }
  
  #compute E(X|W)
  exw = matrix(rep(mux,n),byrow=TRUE,ncol=3)+t(sigxw%*%solve(sigww)%*%t(w-umw))
  protm = exw[,1]; fatm = exw[,2]; energym = exw[,3]
  
  #RC estimation
  rc1 = coxph(Surv(st,event)~protm+fatm+energym+bmi+edu+smoke+factor(centre),data=epic)
  est1 = summary(rc1);
  rc1_res = round(est1$coefficients[1:6,-2],3); 
  
  ###derived method 
  ###protein derived
  fatd1 = fat/prot
  energyd1 = energy/prot
  derivd1 = coxph(Surv(st,event)~prot+fatd1+energyd1+bmi+edu+smoke+factor(centre),data=epic)
  est1 = summary(derivd1);
  derived1_res = est1$coefficients[1:6,-2]
  ##fat derived
  protd2 = prot/fat
  energyd2 = energy/fat
  derivd2 = coxph(Surv(st,event)~protd2+fat+energyd2+bmi+edu+smoke+factor(centre),data=epic)
  est2 = summary(derivd2);
  derived2_res = est2$coefficients[1:6,-2]
  ##energy derived
  protd3 = prot/energy
  fatd3 = fat/energy
  derivd3 = coxph(Surv(st,event)~protd3+fatd3+energy+bmi+edu+smoke+factor(centre),data=epic)
  est3 = summary(derivd3);
  derived3_res = est3$coefficients[1:6,-2]
  ##derived estimation results
  derived_res = rbind(derived1_res[c(1),c(1,2,4)],derived2_res[c(2),c(1,2,4)],
                      derived3_res[c(3),c(1,2,4)],derived2_res[c(4:6),c(1,2,4)])
  
  ###univarite estimation
  cox1 = coxph(Surv(st,event)~prot+bmi+edu+smoke+factor(centre),data=epic)
  lm1 = summary(cox1);
  lm1_res = lm1$coefficients[1:4,-2]
  cox2 = coxph(Surv(st,event)~fat+bmi+edu+smoke+factor(centre),data=epic)
  lm2 = summary(cox2);
  lm2_res = lm2$coefficients[1:4,-2]
  cox3 = coxph(Surv(st,event)~energy+bmi+edu+smoke+factor(centre),data=epic)
  lm3 = summary(cox3);
  lm3_res = lm3$coefficients[1:4,-2]
  univ_res = rbind(lm1_res[1,c(1,2,4)],lm2_res[1,c(1,2,4)],lm3_res[1,c(1,2,4)],
                 lm2_res[2:4,c(1,2,4)])
  res = cbind(naive1_res[,c(1,2,4)],univ_res,derived_res,rc1_res[,c(1,2,4)])
  write.csv(round(res,3),paste0('real_data_ana_sex_',sex,'.csv'))
}

