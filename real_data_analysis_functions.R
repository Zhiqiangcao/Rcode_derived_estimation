# functions used in real data analysis

# The following functions can produce results in Table 5 of manuscript
four_methods_est = function(epic_mask) {
  final_res = NULL
  # analyze male and female separately
  for (sex in 1:2) {
    if (sex == 1)
      epic = epic_mask[epic_mask$sex == 1, ]  #male
    if (sex == 2)
      epic = epic_mask[epic_mask$sex == 2, ]  #female
    
    # sqrt-root transformation, then divided by their standard deviation
    prot = sqrt(epic$qe_prot)
    prot = prot/sd(prot)
    fat = sqrt(epic$qe_fat)
    fat = fat/sd(fat)
    energy = sqrt(epic$qe_energy)
    energy = energy/sd(energy)
    
    ## error-free covariates setup
    bmi = scale(epic$bmi_adj)
    edu = 1 * (epic$l_school >= 3)
    smoke = 1 * (epic$smoke_stat == 3)
    
    ### construct E(x|w) and
    w = cbind(prot, fat, energy)
    sigww = cov(w)  #variance-covariance between variables measured with error
    umw = apply(w, 2, mean)
    n = dim(w)[1]
    
    ########################### naive estimation
    naive1 = coxph(Surv(st, event) ~ prot + fat + energy + bmi + edu + smoke + factor(centre),
                   data = epic)
    est1 = summary(naive1)
    naive1_res = round(est1$coefficients[1:6, -2], 3)
    
    ########################### RC estimation
    if (sex == 1) {
      sigxx = matrix(c(1.05927134, 0.050755983, 0.33722698, 0.05075598, 0.676466609,
                       0.31167934, 0.33722698, 0.31167934, 1.1386985), byrow = T, 3, 3)
      sigxx = round(sigxx, 4)
      siguu = matrix(c(0.8575957, 0.7162641, 0.7517087, 0.7162641, 0.8836357, 0.7924088,
                       0.7517087, 0.7924088, 0.9131195), byrow = T, 3, 3)
      siguu = round(siguu, 4)
      sigxw = matrix(c(0.38838743, 0.021051055, 0.093149182, 0.01860995, 0.280564669,
                       0.086092387, 0.12364605, 0.129269071, 0.31453247), byrow = T, 3, 3)
      sigxw = round(sigxw, 4)
      beta = c(0.3667, 0.4148, 0.2762)
      mux = c(9.356683, 8.654661, 11.937024)
    } else if (sex == 2) {
      sigxx = matrix(c(1.07312719, 0.05953603, 0.32754194, 0.05953603, 1.07539124,
                       0.4410434, 0.32754194, 0.4410434, 1.17795383), byrow = T, 3, 3)
      sigxx = round(sigxx, 4)
      siguu = matrix(c(0.9046743, 0.7202257, 0.7603176, 0.7202257, 0.9308771, 0.8315039,
                       0.7603176, 0.8315039, 0.9254783), byrow = T, 3, 3)
      siguu = round(siguu, 4)
      sigxw = matrix(c(0.31983837, 0.015094115, 0.082384259, 0.01774431, 0.272642958,
                       0.110932463, 0.09762168, 0.11181733, 0.296282223), byrow = T, 3, 3)
      sigxw = round(sigxw, 4)
      beta = c(0.298, 0.2535, 0.2515)
      mux = c(9.245602, 7.510393, 11.167943)
    }
    
    # compute E(X|W)
    exw = matrix(rep(mux, n), byrow = TRUE, ncol = 3) + t(sigxw %*% solve(sigww) %*%
                                                            t(w - umw))
    protm = exw[, 1]
    fatm = exw[, 2]
    energym = exw[, 3]
    
    # RC estimation
    rc1 = coxph(Surv(st, event) ~ protm + fatm + energym + bmi + edu + smoke + factor(centre),
                data = epic)
    est1 = summary(rc1)
    rc1_res = round(est1$coefficients[1:6, -2], 3)
    
    ### derived method protein derived
    fatd1 = fat/prot
    energyd1 = energy/prot
    derivd1 = coxph(Surv(st, event) ~ prot + fatd1 + energyd1 + bmi + edu + smoke + factor(centre),
                    data = epic)
    est1 = summary(derivd1)
    derived1_res = est1$coefficients[1:6, -2]
    ## fat derived
    protd2 = prot/fat
    energyd2 = energy/fat
    derivd2 = coxph(Surv(st, event) ~ protd2 + fat + energyd2 + bmi + edu + smoke + factor(centre),
                    data = epic)
    est2 = summary(derivd2)
    derived2_res = est2$coefficients[1:6, -2]
    ## energy derived
    protd3 = prot/energy
    fatd3 = fat/energy
    derivd3 = coxph(Surv(st, event) ~ protd3 + fatd3 + energy + bmi + edu + smoke + factor(centre),
                    data = epic)
    est3 = summary(derivd3)
    derived3_res = est3$coefficients[1:6, -2]
    ## derived estimation results
    derived_res = rbind(derived1_res[c(1), c(1, 2, 4)], derived2_res[c(2), c(1, 2, 4)],
                        derived3_res[c(3), c(1, 2, 4)], derived2_res[c(4:6), c(1, 2, 4)])
    
    ### univariate estimation
    cox1 = coxph(Surv(st, event) ~ prot + bmi + edu + smoke + factor(centre), data = epic)
    lm1 = summary(cox1)
    lm1_res = lm1$coefficients[1:4, -2]
    cox2 = coxph(Surv(st, event) ~ fat + bmi + edu + smoke + factor(centre), data = epic)
    lm2 = summary(cox2)
    lm2_res = lm2$coefficients[1:4, -2]
    cox3 = coxph(Surv(st, event) ~ energy + bmi + edu + smoke + factor(centre), data = epic)
    lm3 = summary(cox3)
    lm3_res = lm3$coefficients[1:4, -2]
    univ_res = rbind(lm1_res[1, c(1, 2, 4)], lm2_res[1, c(1, 2, 4)], lm3_res[1, c(1,
                                                                                  2, 4)], lm2_res[2:4, c(1, 2, 4)])
    res = cbind(naive1_res[, c(1, 2, 4)], univ_res, derived_res, rc1_res[, c(1, 2, 4)])
    final_res = rbind(final_res, res)
  }
  final_res = round(final_res, 3)
  return(final_res)
}

######The following code can produce results in Table 6 of manuscript##########


# the following procedure is to estimate related results in additive error model 
# (11) of paper, which is realized according to the following reference:
# Zhiqiang Cao and Man Yu Wong. Moment estimation method of parameters in 
# additive measurement error model, Computer Methods and Programs in Biomedicine. 
# 2021, 206:1-9.


### getcfvitac is to estimate unknown parameters for surrogate of vitamin C in
### univariate case, i.e., obtain \sigma^2_{X_c},\sigma^2_{R_c},\sigma^2_{Q_c},
### \sigma^2_{H_c} in Step 1, and correlation coefficient between Q_c and H_c. In
### addition, corresponding correction factor, reliability ratio, var(X_c|Q_c) are also
### calculated Note: Q_c can be FFQ of vegetable+fruit or root vegetable; H_c can be
### 24HR of vegetable+fruit or root vegetable; R_c is plasma vitamin C; X_c is the true
### value of vitamin C

### Input parameters: rho: assumed value of square root of reliability ratio for plasma
### vitamin C (i.e., \rho_c in Step 1) 
### sqchc: sample covariance between Q_c and H_c
### sqcqc: sample variance of Q_c 
### sqcrc: sample covariance between Q_c and R_c 
### shchc: sample variance of H_c 
### shcrc: sample covariance between H_c and R_c 
### srcrc: sample variance of R_c 
### muqc: sample mean of Q_c 
### muhc: sample mean of H_c 
### murc: sample mean of R_c

# output: estimation of unknown parameters for surrogate of vitamin C
getcfvitac = function(rho, sqchc, sqcqc, sqcrc, shchc, shcrc, srcrc, muqc, muhc, murc) {
  sigxcs = rho^2 * srcrc
  sigrcs = srcrc - sigxcs
  sigqcs = sqcqc - sqcrc^2/sigxcs
  sighcs = shchc - shcrc^2/sigxcs
  bqc = sqcrc/sigxcs
  bhc = shcrc/sigxcs
  aqc = muqc - bqc * murc
  ahc = muhc - bhc * murc
  # the following is correlation coefficient between Q_c and H_c
  rho_qchc = (sqchc - bqc * bhc * sigxcs)/sqrt(sigqcs * sighcs)
  # correction factor, reliability ratio, var(X_c|R_c) of plasma vitamin C
  cf_pvc = rho^2
  taus_pvc = sigxcs * (1 - cf_pvc)
  rhos_pvc = cf_pvc
  # correction factor, reliability ratio, var(X_c|Q_c) of vitamin C
  cf_vc = sqcrc/sqcqc
  taus_vc = sigxcs - sqcrc^2/sqcqc
  rhos_vc = sqcrc^2/(sigxcs * sqcqc)
  resu1 = c(aqc, bqc, ahc, bhc, sigxcs, sigrcs, sigqcs, sighcs, rho, rho_qchc)
  names(resu1) = c("aqc", "bqc", "ahc", "bhc", "sigxcs", "sigrcs", "sigqcs", "sighcs", "rho",
                   "rho_qchc")
  resu2 = data.frame(cf = c(cf_pvc, cf_vc), taus = c(taus_pvc, taus_vc), rhos = c(rhos_pvc,
                                                                                  rhos_vc))
  row.names(resu2) = c("plasma vitaminc", "vitaminc")
  resu = list(resu1, resu2)
  return(resu)
}


### getrhoqihi is to find out the range of \rho_{Q_i H_i} such that all estimated
### variances are positive and all estimated correlation coefficients are between -1 and
### 1. Note that ranges of \rho_{X_c X_i}, \rho{Q_c Q_i}, \rho_{Q_c H_i}, \rho_{Q_i
### H_c}, \rho_{H_c H_i}, \sigma^2_{Q_i} and \sigma^2_{H_i} in model (1) of paper are
### calculated (for each value of \rho_{Q_i H_i}). In addition, estimation of correction
### factor, reliability ratio and var(X_i|Q_i) are also calculated Note: Q_i can be FFQ
### of protein, fat or energy; H_i can be 24HR of protein, fat or energy

### Input parameters: 
### sqiqi: sample variance of Q_i 
### shihi: sample variance of H_i 
### sqihi: sample covariance between Q_i and H_i 
### shirc: sample covariance between H_i and R_c
### sqirc: sample covariance between Q_i and R_c 
### sqchc: sample covariance between Q_c and H_c 
### sqcqc: sample variance of Q_c 
### sqcqi: sample covariance between Q_c and Q_i
### sqchi: sample covariance between Q_c and H_i 
### sqcrc: sample covariance between Q_c
### and R_c shchc: sample variance of H_c 
### shcqi: sample covariance between H_c and Q_i
### shchi: sample covariance between H_c and H_i 
### shcrc: sample covariance between H_c and R_c 
### srcrc: sample variance of R_c 
### muqc: sample mean of Q_c 
### muhc: sample mean of H_c 
### murc: sample mean of R_c 
### rho: assumed value of square root of reliability ratio for plasma vitamin C
### bhi: assumed value of \beta_{H_i} in model (1) of paper 
### barh: sample mean of H_i 
### barq: sample mean of Q_i

# output: ranges of \rho_{X_c X_i}, \rho{Q_c Q_i},\rho_{Q_c H_i}, \rho_{Q_i H_c},
# \rho_{H_c H_i}, \sigma^2_{Q_i} and \sigma^2_{H_i} and estimation of correction factor,
# reliability ratio and var(X_i|Q_i)

getrhoqihi = function(sqiqi, shihi, sqihi, shirc, sqirc, sqchc, sqcqc, sqcqi, sqchi, sqcrc,
                      shchc, shcqi, shchi, shcrc, srcrc, muqc, muhc, murc, rho, bhi, barh, barq) {
  # two nonlinear equations
  fun2 <- function(x) {
    f1 <- shihi - (sqihi - rhoqihi * x[1] * x[2]) * (shirc/sqirc) - (x[1])^2
    f2 <- sqiqi - (sqihi - rhoqihi * x[1] * x[2]) * (sqirc/shirc) - (x[2])^2
    return(c(f1 = f1, f2 = f2))
  }
  rhoset = seq(0, 1, by = 0.001)
  nr = length(rhoset)
  sigrhq = matrix(0, nr, 5)
  # default start points
  startx = sqrt(c(shihi/2, sqiqi/2))
  for (i in 1:nr) {
    rhoqihi = rhoset[i]
    resu = multiroot(fun2, startx, maxiter = 500)
    sigrhq[i, ] = round(c(rhoqihi, resu$root, resu$f.root), 6)
  }
  colnames(sigrhq) = c("rhoqihi", "sighi", "sigqi", "froot1", "froot2")
  maxva = sqrt(c(shihi, sqiqi))
  # find out the range of rho_{Q_i H_i}
  i1 = which(sigrhq[, 2] < maxva[1] & sigrhq[, 2] > 0 & abs(sigrhq[, 4]) < 10^(-4))
  i2 = which(sigrhq[, 3] < maxva[2] & sigrhq[, 3] > 0 & abs(sigrhq[, 5]) < 10^(-4))
  minmax = intersect(i1, i2)
  # In this case, cannot obtain target result, need to check data
  if (sum(minmax) == 0) {
    resu = NULL
    return(resu)
  }
  flag = NULL
  for (i in minmax) {
    temp1 = sigrhq[i, 2:3]
    test1 = sum(temp1 > 0)
    test2 = sum(temp1 < maxva)
    temp2 = sigrhq[i, 4:5]
    test3 = sum(temp2 < c(10^(-4), 10^(-4)))
    if (test1 != 2 | test2 != 2 | test3 != 2)
      flag = c(flag, i)
  }
  rorder = setdiff(minmax, flag)
  # In this case, cannot obtain target result
  if (sum(rorder) == 0) {
    resu = NULL
    return(resu)
  }
  resu_temp = sigrhq[rorder, ]
  if (length(rorder) == 1)
    resu_temp = matrix(resu_temp, ncol = 5)
  m = length(rorder)
  rxcxism = rqcqim = rqchim = rqihcm = rhihcm = numeric(m)
  resu0 = getcfvitac(rho, sqchc, sqcqc, sqcrc, shchc, shcrc, srcrc, muqc, muhc, murc)
  resu1 = resu0[[1]]
  bqc = resu1[2]
  bhc = resu1[4]
  sigxcs = resu1[5]
  # \sigma^2_{X_c};
  sigqc = sqrt(resu1[7])
  # \sigma_{Q_c}
  sighc = sqrt(resu1[8])
  # \sigma_{H_c}
  for (i in 1:m) {
    rqihi_i = resu_temp[i, 1]
    sighi_i = resu_temp[i, 2]
    # \sigma_{H_i}
    sigqi_i = resu_temp[i, 3]
    # \sigma_{Q_i}
    rxcxis_i = sqirc * shirc/((sqihi - rqihi_i * sighi_i * sigqi_i) * sigxcs)
    rqcqi_i = (sqcqi - bqc * sqirc)/(sigqc * sigqi_i)
    rqchi_i = (sqchi - bqc * shirc)/(sigqc * sighi_i)
    rqihc_i = (shcqi - bhc * sqirc)/(sighc * sigqi_i)
    rhihc_i = (shchi - bhc * shirc)/(sighc * sighi_i)
    rxcxism[i] = rxcxis_i
    rqcqim[i] = rqcqi_i
    rqchim[i] = rqchi_i
    rqihcm[i] = rqihc_i
    rhihcm[i] = rhihc_i
  }
  i1 = which(sqrt(rxcxism) < 1 & sqrt(rxcxism) > (-1))
  i2 = which(rqcqim < 1 & rqcqim > (-1))
  i3 = which(rqchim < 1 & rqchim > (-1))
  i4 = which(rqihcm < 1 & rqihcm > (-1))
  i5 = which(rhihcm < 1 & rhihcm > (-1))
  rorderf = Reduce(intersect, list(i1, i2, i3, i4, i5))
  if (sum(rorderf) == 0) {
    resu = NULL
    return(resu)
  }
  qihir = resu_temp[rorderf, 1]
  sqhrho1 = sqihi - resu_temp[rorderf, 1] * resu_temp[rorderf, 2] * resu_temp[rorderf, 3]
  temp3 = sqhrho1 * (shirc/sqirc)
  temp4 = sqhrho1 * (sqirc/shirc)
  # range of \beta_{H_i}^2*\sigma^2_{X_i}
  bhixisr = temp3
  # range of \beta_{Q_i}^2*\sigma^2_{X_i}
  bqixisr = temp4
  # range of \hat C.F.
  cf_qir = sqhrho1/(bhi * sqiqi)
  # range of \hat\sigma^2_{X_i}
  sigxisr = temp3/bhi^2
  # range of \hat var(X_i | Q_i)
  tausr = sigxisr - sqhrho1^2/(bhi^2 * sqiqi)
  # range of \hat\rho^2_{Q_i}
  rhosr = sqirc * sqhrho1/(shirc * sqiqi)
  xcxir = sqrt(rxcxism[rorderf])
  # range of \rho_{X_c X_i}
  qcqir = rqcqim[rorderf]
  # range of \rho_{Q_c Q_i}
  qchir = rqchim[rorderf]
  # range of \rho_{Q_c H_i}
  qihcr = rqihcm[rorderf]
  # range of \rho_{Q_i H_c}
  hihcr = rhihcm[rorderf]
  # range of \rho_{H_i H_c}
  sigqisr = resu_temp[rorderf, 3]^2
  # range of \sigma^2_{Q_i}
  sighisr = resu_temp[rorderf, 2]^2
  # range of \sigma^2_{H_i}
  aqi = barq - sqirc/shirc * barh
  aqc = resu1[1]
  bqi = sqirc/shirc * bhi
  nn = length(rorderf)
  resu = data.frame(rqihi = qihir, rxcxi = xcxir, rqcqi = qcqir, rqchi = qchir, rqihc = qihcr,
                    rhihc = hihcr, rsigqis = sigqisr, rsighis = sighisr, rbqixis = bqixisr, rbhixis = bhixisr,
                    rsigxis = sigxisr, rcf = cf_qir, rtaus = tausr, rrhos = rhosr, bqi = rep(bqi, nn),
                    rsqhrho = sqhrho1, aqi = rep(aqi, nn), aqc = rep(aqc, nn), bqc = rep(bqc, nn))
  return(resu)
}

### getrhoqihi_idc is to obtain results from getrhoqihi when one surrogate of true
### vitamin C (i.e., vegetable+fruit or root vegetable) and one nutrient (i.e., protein
### or fat or energy) are specified.

### Input parameters: 
### rho: assumed value of square root of reliability ratio for plasma vitamin C (i.e., \rho_c in Step 1)
### bhi: assumed value of \beta_{H_i} in model (1) of paper 
### idi: one nutrient, which can be protein, fat or energy 
### idc: one surrogate of vitamin C, can be vegetable+fruit
### epic: FFQ and 24HR observed measurements in the validation study of the EPIC-InterAct 

# output: when one surrogate of vitamin C and one nutrient of protein, fat or energy are
# chosen, results are ranges of \rho_{X_c X_i}, \rho{Q_c Q_i},\rho_{Q_c H_i}, \rho_{Q_i
# H_c}, \rho_{H_c H_i}, \sigma^2_{Q_i} and \sigma^2_{H_i} and estimation of correction
# factor, reliability ratio and var(X_i|Q_i)
getrhoqihi_idc = function(rho, bhi, idi, epic) {
  # plasma vitamin c
  Rc = epic$vitC_vitaminC_c
  Rc = Rc/sd(Rc)
  Qc = sqrt(epic$qe_vitc)
  Hc = sqrt(epic$re_vitc)
  Qc = Qc/sd(Qc)
  Hc = Hc/sd(Hc)
  if (idi == 1) {
    Qi = sqrt(epic$qge02 + epic$qge0401)
    Hi = sqrt(epic$rge02 + epic$rge0401)
    Qi = Qi/sd(Qi)
    Hi = Hi/sd(Hi)
  }
  # protein
  if (idi == 2) {
    Qi = sqrt(epic$qe_prot)
    Hi = sqrt(epic$re_prot)
    Qi = Qi/sd(Qi)
    Hi = Hi/sd(Hi)
  }
  # fat
  if (idi == 3) {
    Qi = sqrt(epic$qe_fat)
    Hi = sqrt(epic$re_fat)
    Qi = Qi/sd(Qi)
    Hi = Hi/sd(Hi)
  }
  # energy
  if (idi == 4) {
    Qi = sqrt(epic$qe_energy)
    Hi = sqrt(epic$re_energy)
    Qi = Qi/sd(Qi)
    Hi = Hi/sd(Hi)
  }
  
  # compute summary statistics
  sqcrc = cov(Qc, Rc)
  shcrc = cov(Hc, Rc)
  sqchc = cov(Qc, Hc)
  sqcqc = var(Qc)
  shchc = var(Hc)
  srcrc = var(Rc)
  murc = mean(Rc)
  muqc = mean(Qc)
  muhc = mean(Hc)
  ## one nutrient
  
  sqcqi = cov(Qc, Qi)
  sqchi = cov(Qc, Hi)
  shcqi = cov(Hc, Qi)
  shchi = cov(Hc, Hi)
  sqihi = cov(Qi, Hi)
  sqiqi = var(Qi)
  sqirc = cov(Qi, Rc)
  shirc = cov(Hi, Rc)
  shihi = var(Hi)
  barq1 = mean(Qi)
  barh1 = mean(Hi)
  resu = getrhoqihi(sqiqi, shihi, sqihi, shirc, sqirc, sqchc, sqcqc, sqcqi, sqchi, sqcrc,
                    shchc, shcqi, shchi, shcrc, srcrc, muqc, muhc, murc, rho, bhi, barh1, barq1)
  return(resu)
}

### valbiva_pvitc is to calculate results in C1 of Appendix C when one surrogate of true
### vitamin C and one nutrient of protein, fat or energy are specified

### Input parameters: 
### rho: assumed value of square root of reliability ratio for plasma  vitamin C (i.e., \rho_c in Step 1)
### bhi: assumed value of \beta_{H_i} in model (1) of paper 
### idi: one nutrient, which can be protein, fat or energy 
### idc: one surrogate of vitamin C, can be vegetable+fruit 
### epic: FFQ and 24HR observed measurements in the validation study of the EPIC-InterAct

# output: Results in C1 of Appendix C when one surrogate of true vitamin C and one
# nutrient of protein, fat or energy are specified
valbiva_pvitc = function(rho, bhi, idi, epic) {
  # surrogate of vitamin c is fixed
  resu_idc = getrhoqihi_idc(rho, bhi, idi, epic)
  # range of S_{Q_i H_i}-\rho_{Q_i H_i}*\hat\sigma_{Q_i}*\hat\sigma_{H_i} based on the
  # chosen surrogate of vitamin C
  rsqhrho = resu_idc$rsqhrho
  # find the min and max of ranges of S_{Q_i H_i}-\rho_{Q_i
  # H_i}*\hat\sigma_{Q_i}*\hat\sigma_{H_i}
  minr = min(rsqhrho)
  maxr = max(rsqhrho)
  sqhrho = median(c(minr, maxr))
  # plasma vitamin c
  Rc = epic$vitC_vitaminC_c
  Rc = Rc/sd(Rc)
  ## one nutrient vitamin c
  Qc = sqrt(epic$qe_vitc)
  Hc = sqrt(epic$re_vitc)
  Qc = Qc/sd(Qc)
  Hc = Hc/sd(Hc)
  # fv or root
  if (idi == 1) {
    Qi = sqrt(epic$qge02 + epic$qge0401)
    Hi = sqrt(epic$rge02 + epic$rge0401)
    Qi = Qi/sd(Qi)
    Hi = Hi/sd(Hi)
    umi = mean(Qi/sd(Qi))
    uxi = mean(Hi/sd(Hi))/bhi
  }
  # protein
  if (idi == 2) {
    Qi = sqrt(epic$qe_prot)
    Hi = sqrt(epic$re_prot)
    Qi = Qi/sd(Qi)
    Hi = Hi/sd(Hi)
    umi = mean(Qi/sd(Qi))
    uxi = mean(Hi/sd(Hi))/bhi
  }
  # fat
  if (idi == 3) {
    Qi = sqrt(epic$qe_fat)
    Hi = sqrt(epic$re_fat)
    Qi = Qi/sd(Qi)
    Hi = Hi/sd(Hi)
    umi = mean(Qi/sd(Qi))
    uxi = mean(Hi/sd(Hi))/bhi
  }
  # energy
  if (idi == 4) {
    Qi = sqrt(epic$qe_energy)
    Hi = sqrt(epic$re_energy)
    Qi = Qi/sd(Qi)
    Hi = Hi/sd(Hi)
    umi = mean(Qi/sd(Qi))
    uxi = mean(Hi/sd(Hi))/bhi
  }
  
  umc = mean(Qc/sd(Qc))
  uxc = mean(Rc/sd(Rc))
  um = c(umc, umi)
  ux = c(uxc, uxi)
  
  # compute summary statistics
  sqcrc = cov(Qc, Rc)
  sqcqi = cov(Qc, Qi)
  sqirc = cov(Qi, Rc)
  shirc = cov(Hi, Rc)
  sqcqc = var(Qc)
  srcrc = var(Rc)
  sqiqi = var(Qi)
  
  sigxx = matrix(c(rho^2 * srcrc, shirc/bhi, shirc/bhi, sqhrho * shirc/(bhi^2 * sqirc)),
                 ncol = 2)
  aqi = mean(Qi) - sqirc/shirc * mean(Hi)
  bqi = sqirc/shirc * bhi
  
  
  # variance estimation of error variable for the nutrient
  siguui = sqiqi - sqhrho * sqirc/shirc
  # next compute \hat\Sigma_{WW}, \hat\Sigma_{XW},and \hat\Sigma_{\epsilon \epsilon}
  # of C1 in Appendix C
  
  sigww = matrix(c(sqcqc, sqcqi, sqcqi, sqiqi), ncol = 2)
  sigxw = matrix(c(sqcrc, sqirc, shirc * sqcrc/(bhi * rho^2 * srcrc), sqhrho/bhi), byrow = T,
                 ncol = 2)
  bqc = sqcrc/(rho^2 * srcrc)
  aqc = mean(Qc) - bqc * mean(Rc)
  
  siguu = matrix(c(sqcqc - sqcrc * bqc, sqcqi - sqcrc * sqirc/(rho^2 * srcrc), sqcqi - sqcrc *
                     sqirc/(rho^2 * srcrc), siguui), ncol = 2)
  resu = list(ux, um, sigww, sigxw, sigxx, siguu, aqi, bqi, aqc, bqc, sqhrho, siguui)
  names(resu) = c("mux", "muw", "sigww", "sigxw", "sigxx", "siguu", "aqi", "bqi", "aqc",
                  "bqc", "sqhrho", "siguui")
  return(resu)
}

### getratiof is to compute \rho_{Q_i H_j},\rho_{Q_j H_i}, \rho_{X_i X_j}, \rho_{Q_i
### Q_j} and \rho_{H_i H_j} in Step 3 of Section 2.2, and test whether these five rhos
### are between (-1,1) or not, then, we can narrow results from getrhoqihi_idc based on
### \rho_{Q_i H_i} and \rho_{Q_j H_j} obtained from Step 2. note: where Q_i and Q_j are
### two different nutrients, i.e., Q_i and Q_j can be any two observed values of
### protein, fat and energy

### Input parameters: 
### bhi: assumed value of \beta_{H_i} in model (1) of paper 
### bhj: similar to bhi, i.e., the coefficient of 24HR for another nutrient 
### sigxis: estimated variance of \sigma^2_{X_i} 
### sigxjs: estimated variance of \sigma^2_{X_j} 
### sigqis: estimated variance of \sigma^2_{Q_i} 
### sigqjs: estimated variance of \sigma^2_{Q_j}
### sighis: estimated variance of \sigma^2_{H_i} 
### sighjs: estimated variance of
### \sigma^2_{H_j} sqiqj: sample covariance between Q_i and Q_j 
### sqihj: sample covariance between Q_i and H_j 
### sqjhi: sample covariance between Q_j and H_i 
### shihj: sample covariance between H_i and H_j 
### ci: ci = S_{R_c Q_i} / S_{R_c H_i} 
### cj: cj = S_{R_c Q_j} / S_{R_c H_j} 
### aqh: aqh = \rho_{Q_i H_j} / \rho_{Q_j H_i}, when Q_i, Q_j, H_i and H_j are standardized, 
### we have \rho_{Q_i H_j} = cov(Q_i,H_j), \rho_{Q_j H_i} = cov(Q_j,H_i). 
### Thus, in practice, we used S_{Q_i H_j}, S_{Q_j H_i} to estimate \rho_{Q_i H_j} and \rho_{Q_j H_i}

# output: estimation of \rho_{Q_i H_j},\rho_{Q_j H_i}, \rho_{X_i X_j} \rho_{Q_i Q_j} and
# \rho_{H_i H_j} in Step 3 of Section 2.2, and test whether these five rhos are between
# (-1,1) or not

getratiof = function(bhi, bhj, sigxis, sigxjs, sigqis, sigqjs, sighis, sighjs, sqiqj, sqihj,
                     sqjhi, shihj, ci, cj, aqh) {
  sigxit = sqrt(sigxis)
  sigqit = sqrt(sigqis)
  sighit = sqrt(sighis)
  sigxjt = sqrt(sigxjs)
  sigqjt = sqrt(sigqjs)
  sighjt = sqrt(sighjs)
  rhot = (ci * sqjhi - cj * sqihj)/(ci * sigqjt * sighit - aqh * cj * sigqit * sighjt)
  rqjhimij = rhot
  rqihjmij = aqh * rhot
  tempt = sqihj - aqh * rhot * sigqit * sighjt
  rxixjmij = tempt/(ci * bhi * bhj * sigxit * sigxjt)
  rhihjmij = (ci * shihj - tempt)/(ci * sighit * sighjt)
  rqiqjmij = (sqiqj - cj * tempt)/(sigqit * sigqjt)
  test1 = rqjhimij < 1 & rqjhimij > (-1)
  test2 = rqihjmij < 1 & rqihjmij > (-1)
  test3 = rxixjmij < 1 & rxixjmij > (-1)
  test4 = rhihjmij < 1 & rhihjmij > (-1)
  test5 = rqiqjmij < 1 & rqiqjmij > (-1)
  # test these five rhos are in (-1,1) or not
  test = 1 * (test1 & test2 & test3 & test4 & test5)
  resu = data.frame(xixjrm = rxixjmij, qiqjrm = rqiqjmij, hihjrm = rhihjmij, qihjrm = rqihjmij,
                    qjhirm = rqjhimij, sqihjrhom = tempt, test = test)
  return(resu)
}

### biva is to give the data frame for two chosen nutrients among protein, fat and
### energy

### Input parameters: 
### id1: one nutrient of protein, fat and energy 
### id2: another nutrient
### of protein, fat and energy 
### epic: FFQ and 24HR observed measurements in the validation study of the EPIC-InterAct

# output: which two nutrients are chosen
biva = function(id1, id2, epic) {
  if (id1 == 1) {
    # fv+ fv+prot
    if (id2 == 2) {
      Qi = sqrt(epic$qge02 + epic$qge0401)
      Qj = sqrt(epic$qe_prot)
      Hi = sqrt(epic$rge02 + epic$rge0401)
      Hj = sqrt(epic$re_prot)
      Qi = Qi/sd(Qi)
      Hi = Hi/sd(Hi)
      Qj = Qj/sd(Qj)
      Hj = Hj/sd(Hj)
    }
    if (id2 == 3) {
      Qi = sqrt(epic$qge02 + epic$qge0401)
      Qj = sqrt(epic$qe_fat)
      Hi = sqrt(epic$rge02 + epic$rge0401)
      Hj = sqrt(epic$re_fat)
      Qi = Qi/sd(Qi)
      Hi = Hi/sd(Hi)
      Qj = Qj/sd(Qj)
      Hj = Hj/sd(Hj)
    }
    if (id2 == 4) {
      Qi = sqrt(epic$qge02 + epic$qge0401)
      Qj = sqrt(epic$qe_energy)
      Hi = sqrt(epic$rge02 + epic$rge0401)
      Hj = sqrt(epic$re_energy)
      Qi = Qi/sd(Qi)
      Hi = Hi/sd(Hi)
      Qj = Qj/sd(Qj)
      Hj = Hj/sd(Hj)
    }
  }
  
  if (id1 == 2) {
    ## prot+ prot+fat
    if (id2 == 3) {
      Qi = sqrt(epic$qe_prot)
      Qj = sqrt(epic$qe_fat)
      Hi = sqrt(epic$re_prot)
      Hj = sqrt(epic$re_fat)
      Qi = Qi/sd(Qi)
      Hi = Hi/sd(Hi)
      Qj = Qj/sd(Qj)
      Hj = Hj/sd(Hj)
    }
    # prot+energy
    if (id2 == 4) {
      Qi = sqrt(epic$qe_prot)
      Qj = sqrt(epic$qe_energy)
      Hi = sqrt(epic$re_prot)
      Hj = sqrt(epic$re_energy)
      Qi = Qi/sd(Qi)
      Hi = Hi/sd(Hi)
      Qj = Qj/sd(Qj)
      Hj = Hj/sd(Hj)
    }
  }
  if (id1 == 3) {
    ## fat+ fat+energy
    if (id2 == 4) {
      Qi = sqrt(epic$qe_fat)
      Qj = sqrt(epic$qe_energy)
      Hi = sqrt(epic$re_fat)
      Hj = sqrt(epic$re_energy)
      Qi = Qi/sd(Qi)
      Hi = Hi/sd(Hi)
      Qj = Qj/sd(Qj)
      Hj = Hj/sd(Hj)
    }
  }
  
  resu = data.frame(Qi, Hi, Qj, Hj)
  return(resu)
}

### valbiva_idc is to calculate results from getrhoqihi_idc for the same surrogate of
### true vitamin C and two different nutrients of protein, fat and energy

### Input parameters: 
### rho: assumed value of square root of reliability ratio for plasma vitamin C 
### bhi: assumed value of \beta_{H_i} in model (1) of paper 
### bhj: similar to bhi, i.e., the coefficient of 24Hr for another nutrient 
### id1: one nutrient of protein, fat or energy 
### id2: another nutrient of protein, fat and energy 
### idc: one surrogate of vitamin C, can be vegetable+fruit
### epic: FFQ and 24HR observed measurements in the validation study of the EPIC-InterAct

# output: return of results from getrhoqihi_idc for the same surrogate of vitamin C and
# two different nutrients of protien, fat and energy, observed values of these two
# different nutrients are also returned

valbiva_idc = function(rho, bhi, bhj, id1, id2, epic) {
  # check which combination of nutrients is used
  qhresu = biva(id1, id2, epic)
  resui = getrhoqihi_idc(rho, bhi, idi = id1, epic)
  resuj = getrhoqihi_idc(rho, bhj, idi = id2, epic)
  resu = list(resui, resuj, qhresu)
  return(resu)
}

### getratioff is to compute values of S_{Q_i H_j}-aqh*\hat\rho_{Q_j
### H_i}*\hat\sigma_{Q_i}*\hat\sigma_{H_j} for all combinations with one from
### \sigma_{Q_i}, \sigma_{H_i}, \sigma_{X_i} and one from \sigma_{Q_j}, \sigma_{H_j},
### \sigma_{X_j}, then record whether each value is reasonable or not (i.e., if all five
### rhos in getratiof are between (-1,1), then the corresponding S_{Q_i
### H_j}-aqh*\hat\rho_{Q_j H_i}*\hat\sigma_{Q_i}*\hat\sigma_{H_j} is said to be
### reasonable)

### Input parameters: 
### bhi: assumed value of \beta_{H_i} in model (1) of paper 
### bhj: similar to bhi, i.e., the coefficient of 24HR for another nutrient 
### resui: the first return object of subprogram valbiva_idc 
### resuj: the second return object of subprogram valbiva_idc 
### index1: index of common range of S_{Q_i H_i}-\rho_{Q_i H_i}\hat\sigma_{Q_i}\hat\sigma_{H_i}
### when one nutrient is used 
### index2: index of common range of S_{Q_j H_j}-\rho_{Q_j H_j}\hat\sigma_{Q_j}\hat\sigma_{H_j} when
### another nutrient is used 
### sqiqj: sample covariance between Q_i and Q_j 
### sqihj: sample covariance between Q_i and H_j 
### sqjhi: sample covariance between Q_j and H_i 
### shihi: sample covariance between H_i and H_j 
### ci: ci = S_{R_c Q_i} / S_{R_c H_i} 
### cj: cj = S_{R_c Q_j} / S_{R_c H_j} 
### aqh: aqh = \rho_{Q_i H_j} / \rho_{Q_j H_i}, when Q_i, Q_j,
### H_i and H_j are standardized, we have \rho_{Q_i H_j} = cov(Q_i,H_j), 
### \rho_{Q_j H_i} = cov(Q_j,H_i). Thus, in practice, we use S_{Q_i H_j}, 
### S_{Q_j H_i} to estimate \rho_{Q_i H_j} and \rho_{Q_j H_i}

# output: values of S_{Q_i H_j}-aqh*\hat\rho_{Q_j
# H_i}*\hat\sigma_{Q_i}*\hat\sigma_{H_j}, as well as an indicator whether these values
# are reasonable or not

getratioff = function(bhi, bhj, resui, resuj, index1, index2, sqiqj, sqihj, sqjhi, shihj,
                      ci, cj, aqh) {
  resui = resui[index1, ]
  resuj = resuj[index2, ]
  sigxis = resui[, 11]
  sigqis = resui[, 7]
  sighis = resui[, 8]
  sigxjs = resuj[, 11]
  sigqjs = resuj[, 7]
  sighjs = resuj[, 8]
  nn1 = length(sigxis)
  nn2 = length(sigxjs)
  rhoij = NULL
  test = NULL
  for (i in 1:nn1) {
    sigxisi = sigxis[i]
    sigqisi = sigqis[i]
    sighisi = sighis[i]
    tempresu = getratiof(bhi, bhj, sigxisi, sigxjs, sigqisi, sigqjs, sighisi, sighjs,
                         sqiqj, sqihj, sqjhi, shihj, ci, cj, aqh)
    rhoij = c(rhoij, tempresu$sqihjrhom)
    test = c(test, tempresu$test)
  }
  resu = data.frame(rhoij = rhoij, test = test)
  return(resu)
}


### valbiva_onut is to find corresponding estimated results in C2 of Appendix C for
### surrogates for any two true nutrients of protein, fat and energy

### Input parameters: 
### rho: assumed value of square root of reliability ratio for plasma vitamin C (i.e., \rho_c in Step 1) 
### bhi: assumed value of \beta_{H_i} in model (1) of paper 
### bhj: similar to bhi, i.e., the coefficient of 24HR for another nutrient 
### id1: one nutrient of protein, fat or energy 
### id2: another nutrient of protein, fat or energy 
### idc: one surrogate of vitamin C, can be vegetable+fruit
### epic: FFQ and 24HR observed measurements in the validation study of the EPIC-InterAct
### len: default value is 101, which is for the length of final common range of
### S_{Q_i H_j}-aqh*\hat \rho_{Q_j H_i}*\hat\sigma_{Q_i}*\hat\sigma_{H_j}

# output: corresponding estimated results in C2 of Appendix when two different nutrients
# of protein, fat and energy are chosen
valbiva_onut = function(rho, bhi, bhj, id1, id2, epic, len = 101) {
  # vegetable+fruit or root vegetable as surrogate of vitamin c
  resu_valbiva_idc3 = valbiva_idc(rho, bhi, bhj, id1, id2, epic)
  
  # next, find the common range of S_{Q_i H_i}-\rho_{Q_i
  # H_i}\hat\sigma_{Q_i}\hat\sigma_{H_i} based on idc=1 and idc=2
  rsqhrhoi3 = resu_valbiva_idc3[[1]]$rsqhrho
  minrii = min(rsqhrhoi3)
  maxrii = max(rsqhrhoi3)
  # find the common range of S_{Q_j H_j}-\rho_{Q_j
  # H_j}\hat\sigma_{Q_j}\hat\sigma_{H_j} based on idc=1 and idc=2
  rsqhrhoj3 = resu_valbiva_idc3[[2]]$rsqhrho
  minrjj = min(rsqhrhoj3)
  maxrjj = max(rsqhrhoj3)
  # next, find the common range for S_{Q_i H_j}-aqh*\hat\rho_{Q_j
  # H_i}*\hat\sigma_{Q_i}*\hat\sigma_{H_j}
  qhresu = resu_valbiva_idc3[[3]]  #two different nutrients
  Rc = epic$vitC_vitaminC_c
  Qi = qhresu$Qi
  Qj = qhresu$Qj
  Hi = qhresu$Hi
  Hj = qhresu$Hj
  # corresponding summary statistics
  sqiqi = var(Qi)
  sqjqj = var(Qj)
  shihi = var(Hi)
  shjhj = var(Hj)
  sqirc = cov(Qi, Rc)
  shirc = cov(Hi, Rc)
  sqjqj = var(Qj)
  sqjrc = cov(Qj, Rc)
  shjrc = cov(Hj, Rc)
  sqiqj = cov(Qi, Qj)
  shihj = cov(Hi, Hj)
  sqihj = cov(Qi, Hj)
  sqjhi = cov(Qj, Hi)
  ci = sqirc/shirc
  cj = sqjrc/shjrc
  bqi = (sqirc/shirc) * bhi
  aqi = mean(Qi) - (sqirc/shirc) * mean(Hi)
  bqj = (sqjrc/shjrc) * bhj
  aqj = mean(Qj) - (sqjrc/shjrc) * mean(Hj)
  emratio = sqihj/sqjhi
  # aqh = cov(Q_i,H_j)/cov(Q_j,H_i) when Q_i, H_j, Q_j and H_i are standardized Here,
  # we use sqihj=S_{Q_i H_j}, sqjhi=S_{Q_j H_i} to estimate \rho_{Q_i,H_j} and
  # \rho_{Q_j,H_i} find the common range of rhoqihj
  indexi3 = which(rsqhrhoi3 >= minrii & rsqhrhoi3 <= maxrii)
  indexj3 = which(rsqhrhoj3 >= minrjj & rsqhrhoj3 <= maxrjj)
  # based on common range of S_{Q_i H_i}-\rho_{Q_i
  # H_i}\hat\sigma_{Q_i}\hat\sigma_{H_i} to find valid range of S_{Q_i
  # H_j}-aqh*\hat\rho_{Q_j H_i}*\hat\sigma_{Q_i}*\hat\sigma_{H_j}
  rhoresu3 = getratioff(bhi, bhj, resu_valbiva_idc3[[1]], resu_valbiva_idc3[[2]], indexi3,
                        indexj3, sqiqj, sqihj, sqjhi, shihj, ci, cj, aqh = emratio)
  rsqhrhoij3 = rhoresu3$rhoij[rhoresu3$test == 1]
  
  minrij = min(rsqhrhoij3)
  maxrij = max(rsqhrhoij3)
  # result of rsqhrhoi,rsqhrhoj are obtained when choice is determined the following
  # is the fourth assumption
  sqhrhoi = median(c(minrii, maxrii))
  sqhrhoj = median(c(minrjj, maxrjj))
  # variance and covariance between Qi and Qj
  sigxis = sqhrhoi/(ci * bhi^2)
  sigxjs = sqhrhoj/(cj * bhj^2)
  sigqis = sqiqi - ci * sqhrhoi
  sigqjs = sqjqj - cj * sqhrhoj
  # sqihjrho=temsure$sqihjrhom, 2*1 matrix, assumed \alpha_{h_i}=\alpha_{h_j}=0
  ux = matrix(c(mean(Hi)/bhi, mean(Hj)/bhj), ncol = 1)
  # len = length(sqihjrhoset)
  sigww = matrix(c(sqiqi, sqiqj, sqiqj, sqjqj), ncol = 2)
  # quantile of the final range of S_{Q_i H_j}-aqh*\hat\rho_{Q_j
  # H_i}*\hat\sigma_{Q_i}*\hat\sigma_{H_j} here, we set the length of final range is
  # 101. and record which values of S_{Q_i H_j}-aqh*\hat\rho_{Q_j
  # H_i}*\hat\sigma_{Q_i}*\hat\sigma_{H_j} can produce valid matrix of \Sigma_{XX} and
  # \Sigma_{UU} (i.e., \Sigma_{\epsilon \epsilon} in C2 of Appendix C)
  sqihjrhoset = quantile(c(minrij, maxrij), probs = seq(0, 1, length = len))
  sigxxm = siguum = sigxwm = array(0, dim = c(2, 2, len))
  rhoxu = matrix(0, len, 3)
  for (i in 1:len) {
    sqihjrho = sqihjrhoset[i]
    covxij = sqihjrho/(ci * bhi * bhj)
    covuij = sqiqj - cj * sqihjrho
    sigxw = matrix(c(sqhrhoi/bhi, cj * sqihjrho/(ci * bhi), sqihjrho/bhj, sqhrhoj/bhj),
                   byrow = T, ncol = 2)
    sigxx = matrix(c(sigxis, covxij, covxij, sigxjs), ncol = 2)
    siguu = matrix(c(sigqis, covuij, covuij, sigqjs), ncol = 2)
    # if rhoxu[i,3]=1, then estimated \sigma_{XX} and \sigma_{UU} are valid
    if (det(sigxx) > 0 & det(siguu) > 0)
      rhoxu[i, 3] = 1
    sigxxm[, , i] = sigxx
    siguum[, , i] = siguu
    sigxwm[, , i] = sigxw
    rhoxu[i, 1] = sigxx[1, 2]/sqrt(prod(diag(sigxx)))
    rhoxu[i, 2] = siguu[1, 2]/sqrt(prod(diag(siguu)))
  }
  rhoxu = cbind(rhoxu, sqihjrhoset)
  colnames(rhoxu) = c("rhox1x2", "rhou1u2", "True or False", "rhoqihj")
  resu = list(ux, sigww, sigxwm, sigxxm, siguum, aqi, bqi, aqj, bqj, minrij, maxrij, rhoxu)
  names(resu) = c("mux", "sigww", "sigxw", "sigxx", "siguu", "aqi", "bqi", "aqj", "bqj",
                  "minrij", "maxrij", "rhoxu")
  return(resu)
}

### valid_resu is to compute the final estimation results similar to the example in
### Section 3.2 of paper, i.e., four covariates measured with error, W =
### (vegetable+fruit, protein, fat, energy)^T or W = (root vegetable, protein, fat,
### energy)^T note: for estimating \Sigma_{XX},\Sigma_{WW},\Sigma_{XW},\Sigma_{\epsilon
### \epsilon}, this program first calculates all bivariate estimation results (C1 and C2
### of Appendix C), then construct estimates of the four 4*4 matrices of
### \Sigma_{XX},\Sigma_{WW},\Sigma_{XW},\Sigma_{\epsilon \epsilon} according to their
### bivariate estimation results.


### Input parameters: 
### epic: FFQ and 24HR observed measurements in the validation study  of the EPIC-InterAct
### sex: sex=1 is for male and sex=2 is for female 
### rho: assumed value of square root of reliability ratio for plasma vitamin C (i.e., \rho_c in Step 1) 
### bhi: assumed value of \beta_{H_i} in model (1) of paper 
### bhj: similar to bhi, i.e., the coefficient of 24HR for another nutrient

# output: estimation results of four covariates measured with error, i.e., estimates of
# \Sigma_{XX}, \Sigma_{XW}, \Sigma_{\epsilon \epsilon} similar to the example in Section
# 3.2 of paper

valid_resu = function(epic, sex, rho, bhi, bhj) {
  # sex takes 1 or 2
  epic = epic[epic$sex == sex, ]
  len = 101
  # results for surrogates of true vitamin C and protein in C1 of Appendix C
  resu = valbiva_pvitc(rho, bhi, idi = 1, epic)
  sigxx_vcc = resu$sigxx
  siguu_vcc = resu$siguu
  sigxw_vcc = resu$sigxw
  bqi_vcc = resu$bqi
  bqc_vcc = resu$bqc
  mux_vcc = resu$mux
  muw_vcc = resu$muw
  
  # results for surrogates of true vitamin C and protein in C1 of Appendix C
  resu = valbiva_pvitc(rho, bhi, idi = 2, epic)
  sigxx_vcp = resu$sigxx
  siguu_vcp = resu$siguu
  sigxw_vcp = resu$sigxw
  bqi_vcp = resu$bqi
  bqc_vcp = resu$bqc
  mux_vcp = resu$mux
  muw_vcp = resu$muw
  
  # results for surrogates of true vitamin C and fat in C1 of Appendix C
  resu = valbiva_pvitc(rho, bhi, idi = 3, epic)
  sigxx_vcf = resu$sigxx
  siguu_vcf = resu$siguu
  sigxw_vcf = resu$sigxw
  bqi_vcf = resu$bqi
  bqc_vcf = resu$bqc
  mu_vcf = resu$mux
  mux_vcf = resu$mux
  muw_vcf = resu$muw
  # results for surrogates of true vitamin C and energy in C1 of Appendix C
  resu = valbiva_pvitc(rho, bhi, idi = 4, epic)
  sigxx_vce = resu$sigxx
  siguu_vce = resu$siguu
  sigxw_vce = resu$sigxw
  bqi_vce = resu$bqi
  bqc_vce = resu$bqc
  mux_vce = resu$mux
  muw_vce = resu$muw
  
  # results for surrogates of protein and fat in in C2 of Appendix C
  resu = valbiva_onut(rho, bhi, bhj, id1 = 1, id2 = 2, epic, len)
  rhoxu = resu$rhoxu
  index = rhoxu[, 3]
  sigxx = resu$sigxx[, , index == 1]
  # keep those valid estimation results
  siguu = resu$siguu[, , index == 1]
  sigxw = resu$sigxw[, , index == 1]
  sigxxn = siguun = sigxwn = matrix(0, 2, 2)
  len1 = sum(index)
  for (i in 1:len1) {
    sigxxn = sigxxn + sigxx[, , i]
    siguun = siguun + siguu[, , i]
    sigxwn = sigxwn + sigxw[, , i]
  }
  sigx12 = sigxxn/len1
  # use the mean of those valid estimated results
  sigu12 = siguun/len1
  sigxw12 = sigxwn/len1
  sigww12 = resu$sigww
  
  resu = valbiva_onut(rho, bhi, bhj, id1 = 1, id2 = 3, epic, len)
  rhoxu = resu$rhoxu
  index = rhoxu[, 3]
  sigxx = resu$sigxx[, , index == 1]
  siguu = resu$siguu[, , index == 1]
  sigxw = resu$sigxw[, , index == 1]
  sigxxn = siguun = sigxwn = matrix(0, 2, 2)
  len1 = sum(index)
  for (i in 1:len1) {
    sigxxn = sigxxn + sigxx[, , i]
    siguun = siguun + siguu[, , i]
    sigxwn = sigxwn + sigxw[, , i]
  }
  sigx13 = sigxxn/len1
  sigu13 = siguun/len1
  sigxw13 = sigxwn/len1
  sigww13 = resu$sigww
  
  resu = valbiva_onut(rho, bhi, bhj, id1 = 1, id2 = 4, epic, len)
  rhoxu = resu$rhoxu
  index = rhoxu[, 3]
  sigxx = resu$sigxx[, , index == 1]
  siguu = resu$siguu[, , index == 1]
  sigxw = resu$sigxw[, , index == 1]
  sigxxn = siguun = sigxwn = matrix(0, 2, 2)
  
  len1 = sum(index)
  for (i in 1:len1) {
    sigxxn = sigxxn + sigxx[, , i]
    siguun = siguun + siguu[, , i]
    sigxwn = sigxwn + sigxw[, , i]
  }
  sigx14 = sigxxn/len1
  sigu14 = siguun/len1
  sigxw14 = sigxwn/len1
  sigww14 = resu$sigww
  
  # results for surrogates of fat and energy in C2 of Appendix C
  resu = valbiva_onut(rho, bhi, bhj, id1 = 2, id2 = 3, epic, len)
  rhoxu = resu$rhoxu
  index = rhoxu[, 3]
  sigxx = resu$sigxx[, , index == 1]
  siguu = resu$siguu[, , index == 1]
  sigxw = resu$sigxw[, , index == 1]
  sigxxn = siguun = sigxwn = matrix(0, 2, 2)
  
  len1 = sum(index)
  for (i in 1:len1) {
    sigxxn = sigxxn + sigxx[, , i]
    siguun = siguun + siguu[, , i]
    sigxwn = sigxwn + sigxw[, , i]
  }
  sigx23 = sigxxn/len1
  sigu23 = siguun/len1
  sigxw23 = sigxwn/len1
  sigww23 = resu$sigww
  # results for surrogates of fat and energy in C2 of Appendix C
  resu = valbiva_onut(rho, bhi, bhj, id1 = 2, id2 = 4, epic, len)
  rhoxu = resu$rhoxu
  index = rhoxu[, 3]
  sigxx = resu$sigxx[, , index == 1]
  siguu = resu$siguu[, , index == 1]
  sigxw = resu$sigxw[, , index == 1]
  sigxxn = siguun = sigxwn = matrix(0, 2, 2)
  
  len1 = sum(index)
  for (i in 1:len1) {
    sigxxn = sigxxn + sigxx[, , i]
    siguun = siguun + siguu[, , i]
    sigxwn = sigxwn + sigxw[, , i]
  }
  sigx24 = sigxxn/len1
  sigu24 = siguun/len1
  sigxw24 = sigxwn/len1
  sigww24 = resu$sigww
  
  # results for surrogates of fat and energy in C2 of Appendix C
  resu = valbiva_onut(rho, bhi, bhj, id1 = 3, id2 = 4, epic, len)
  rhoxu = resu$rhoxu
  index = rhoxu[, 3]
  sigxx = resu$sigxx[, , index == 1]
  siguu = resu$siguu[, , index == 1]
  sigxw = resu$sigxw[, , index == 1]
  sigxxn = siguun = sigxwn = matrix(0, 2, 2)
  len1 = sum(index)
  for (i in 1:len1) {
    sigxxn = sigxxn + sigxx[, , i]
    siguun = siguun + siguu[, , i]
    sigxwn = sigxwn + sigxw[, , i]
  }
  sigx34 = sigxxn/len1
  sigu34 = siguun/len1
  sigxw34 = sigxwn/len1
  sigww34 = resu$sigww
  
  # four covariates
  beta_hat = c(bqi_vcp, bqi_vcf, bqi_vce)
  mux_hat = c(mux_vcp[2], mux_vcf[2], mux_vce[2])
  
  # the following three matrices are for protein, fat and energy
  sigxx4 = matrix(c(sigx12[1, 1], sigx12[1, 2], sigx13[1, 2], sigx14[1, 2], sigx12[1, 2],
                    sigx12[2, 2], sigx23[1, 2], sigx24[1, 2], sigx13[1, 2], sigx23[1, 2], sigx23[2, 2],
                    sigx34[1, 2], sigx14[1, 2], sigx24[1, 2], sigx34[1, 2], sigx34[2, 2]), byrow = T,
                  4, 4)
  siguu4 = matrix(c(sigu12[1, 1], sigu12[1, 2], sigu13[1, 2], sigu14[1, 2], sigu12[1, 2],
                    sigu12[2, 2], sigu23[1, 2], sigu24[1, 2], sigu13[1, 2], sigu23[1, 2], sigu23[2, 2],
                    sigu34[1, 2], sigu14[1, 2], sigu24[1, 2], sigu34[1, 2], sigu34[2, 2]), byrow = T,
                  4, 4)
  
  resu = list(mux_hat, beta_hat, sigxx4, siguu4)
  names(resu) = c("mut", "b", "Sigma_T", "Sigma_e")
  return(resu)
}

# compute correlation coefficients based on a variance-covariance matrix
cal_corr_matrix = function(x) {
  # note: here x is a 3*3 matrix
  x[1, 2] = x[2, 1] = x[1, 2]/sqrt(x[1, 1] * x[2, 2])
  x[1, 3] = x[3, 1] = x[1, 3]/sqrt(x[1, 1] * x[3, 3])
  x[2, 3] = x[3, 2] = x[2, 3]/sqrt(x[2, 2] * x[3, 3])
  return(x)
}
