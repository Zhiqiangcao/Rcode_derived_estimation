### This R code shows how the main study data and validation study data from EPIC-InterAct Study is masked

epic = read.csv("./mainstudyf_2.csv", header = TRUE)
# mainstudyf_2.csv is real data of main study in EPIC-InterAct Study, which is
# confidential
epic1 = epic[, c(4, 13:16, 47, 46, 36:37, 39)]
n = dim(epic1)[1]

# generate noise
set.seed(123456)
bmi_noise = rnorm(n, sd = 0.5)
prot_noise = rnorm(n, sd = 5)
fat_noise = rnorm(n, sd = 5)
energy_noise = rnorm(n, sd = 50)

# added noise for those continuous variables used in paper
epic1$bmi_adj = epic1$bmi_adj + bmi_noise
epic1$qe_prot = epic1$qe_prot + prot_noise
epic1$qe_fat = epic1$qe_fat + fat_noise
epic1$qe_energy = epic1$qe_energy + energy_noise

# delete those observations with less than zero values and then delete random 200
# samples
epic2 = epic1[-2237, ]
set.seed(123456)
ni = sample(1:15256, 200)
epic3 = epic2[-ni, ]
dim(epic3)
# [1] 15056 10 
# note: the sample size of the masked epic data is 15056, which is
# different from that in paper
epic_mask = epic3

save(epic_mask, file = "./epic_mask.RData")

########### Validation study
validation = read.csv("./validationf.csv", header = TRUE)
#validationf.csv is real data of validation study dataset in EPIC-InterAct Study, which is confidential
validation1 = validation[, c(13, 17, 27, 36:37,39, 41, 42, 43, 53, 62:63,65,67)] #keep variables used in validation study
n = dim(validation1)[1]

# generate noise
set.seed(123456)
plasma_vitac_noise = runif(n)
qeprot_noise = rnorm(n, sd = 3)
reprot_noise = rnorm(n, sd = 2.5)
qefat_noise = rnorm(n, sd = 3)
refat_noise = rnorm(n, sd = 2.5)
qeenergy_noise = rnorm(n, sd = 3)
reenergy_noise = rnorm(n, sd = 2.5)
qge02_noise = rnorm(n, sd = 3)
rge02_noise = rnorm(n, sd = 3.5)
qge0401_noise = rnorm(n, sd = 3)
rge0401_noise = rnorm(n, sd = 3.5)
qevitc_noise = rnorm(n, sd = 3)
revitc_noise = rnorm(n, sd = 2.5)

# added noise
validation1$vitC_vitaminC_c = validation1$vitC_vitaminC_c + plasma_vitac_noise
validation1$qe_prot = validation1$qe_prot + qeprot_noise
validation1$re_prot = validation1$re_prot + reprot_noise
validation1$qe_fat = validation1$qe_fat + qefat_noise
validation1$re_fat = validation1$re_fat + refat_noise
validation1$qe_energy = validation1$qe_energy + qeenergy_noise
validation1$re_energy = validation1$re_energy + reenergy_noise
validation1$qge02 = validation1$qge02 + qge02_noise
validation1$rge02 = validation1$rge02 + rge02_noise
validation1$qge0401 = validation1$qge0401 + qge0401_noise
validation1$rge0401 = validation1$rge0401 + rge0401_noise
validation1$qe_vitc = validation1$qe_vitc + qevitc_noise
validation1$re_vitc = validation1$re_vitc + revitc_noise

# delete observations with unreasonable values (i.e., less than zero values)
index1 = validation1$qge02 < 0
index2 = validation1$rge02 < 0
index3 = validation1$qge0401 < 0
index4 = validation1$rge0401 < 0
index5 = validation1$re_vitc <0 
index = index1 | index2 | index3 | index4 | index5 
validation2 = validation1[!index, ]
dim(validation2)
#870 14

# note: the sample size of the masked validation study is 870, which is
# different from that in paper

validation_mask = validation2
save(validation_mask, file = "./validation_mask.RData")
