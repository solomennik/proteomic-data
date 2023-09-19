########################################################################################################################################################################################################################
#Load Packages
########################################################################################################################################################################################################################
library(truncnorm)
library(lattice)
library(nlme)
library(MCMCpack)
library(combinat)

########################################################################################################################################################################################################################
#Load Data
########################################################################################################################################################################################################################
data=read.csv("Training_data.csv")
data_validation=read.csv("Validation_data.csv")

########################################################################################################################################################################################################################
#Fit linear mixed models in training data controls to get estimates of mu_theta_k, sigma_k_2 and sigma_theta_k_2 (pg 14 Supplementary Materials)
########################################################################################################################################################################################################################
control_data=subset(data,data$D==0)

hist(control_data$Y1)
hist(control_data$Y2)
hist(control_data$Y3)

#Fit linear mixed models with random intercept
model_1=lme(Y1~1, random=~1| ID ,data=control_data)
summary(model_1)
model_2=lme(Y2~1, random=~1| ID ,data=control_data)
summary(model_2)
model_3=lme(Y3~1, random=~1| ID ,data=control_data)
summary(model_3)

########################################################################################################################################################################################################################
#Implement PEB algorithm in validation dataset
########################################################################################################################################################################################################################

#Biomarker 1
mu_theta_1=model_1$coefficients$fixed
sigma_1_2=as.numeric(VarCorr(model_1)[2,1])
sigma_theta_1_2=as.numeric(VarCorr(model_1)[1,1])
B_1=sigma_theta_1_2/(sigma_1_2+sigma_theta_1_2)
Z1=(data_validation$Y1-mu_theta_1)/sqrt(sigma_1_2+sigma_theta_1_2)
B_n_1=sigma_theta_1_2/((sigma_1_2/(data_validation$obs_number-1))+sigma_theta_1_2)

#Biomarker 2
mu_theta_2=model_2$coefficients$fixed
sigma_2_2=as.numeric(VarCorr(model_2)[2,1])
sigma_theta_2_2=as.numeric(VarCorr(model_2)[1,1])
B_2=sigma_theta_2_2/(sigma_2_2+sigma_theta_2_2)
Z2=(data_validation$Y2-mu_theta_2)/sqrt(sigma_2_2+sigma_theta_2_2)
B_n_2=sigma_theta_2_2/((sigma_2_2/(data_validation$obs_number-1))+sigma_theta_2_2)

#Biomarker 3
mu_theta_3=model_3$coefficients$fixed
sigma_3_2=as.numeric(VarCorr(model_3)[2,1])
sigma_theta_3_2=as.numeric(VarCorr(model_3)[1,1])
B_3=sigma_theta_3_2/(sigma_3_2+sigma_theta_3_2)
Z3=(data_validation$Y3-mu_theta_3)/sqrt(sigma_3_2+sigma_theta_3_2)
B_n_3=sigma_theta_3_2/((sigma_3_2/(data_validation$obs_number-1))+sigma_theta_3_2)

#Get sample size 
N=length(unique(data_validation$ID)) #total number of patients
n_0=length(unique(data_validation$ID[(data_validation$D==0)])) #number of control patients
subj_ID=unique(data_validation$ID)

Z1_mean=rep(NA,length(data_validation$ID))
Z2_mean=rep(NA,length(data_validation$ID))
Z3_mean=rep(NA,length(data_validation$ID))
for(i in 1:length(subj_ID))
{
  #biomarker 1
  Z_subject=Z1[data_validation$ID==subj_ID[i]]
  temp=array(Z_subject,c(length(Z_subject),length(Z_subject)))
  temp[lower.tri(temp,diag=TRUE)]=NA
  Z1_mean[data_validation$ID==subj_ID[i]]=apply(temp,2,mean,na.rm=TRUE)
  
  #biomarker 2
  Z_subject=Z2[data_validation$ID==subj_ID[i]]
  temp=array(Z_subject,c(length(Z_subject),length(Z_subject)))
  temp[lower.tri(temp,diag=TRUE)]=NA
  Z2_mean[data_validation$ID==subj_ID[i]]=apply(temp,2,mean,na.rm=TRUE)
  
  #biomarker 3
  Z_subject=Z3[data_validation$ID==subj_ID[i]]
  temp=array(Z_subject,c(length(Z_subject),length(Z_subject)))
  temp[lower.tri(temp,diag=TRUE)]=NA
  Z3_mean[data_validation$ID==subj_ID[i]]=apply(temp,2,mean,na.rm=TRUE)
}
Z1_mean[(data_validation$obs_number==1)]=0
Z2_mean[(data_validation$obs_number==1)]=0
Z3_mean[(data_validation$obs_number==1)]=0

data_validation$Z1_PEB=(Z1-Z1_mean*B_n_1)/sqrt(1-B_1*B_n_1)
data_validation$Z2_PEB=(Z2-Z2_mean*B_n_2)/sqrt(1-B_2*B_n_2)
data_validation$Z3_PEB=(Z3-Z3_mean*B_n_3)/sqrt(1-B_3*B_n_3)

########################################################################################################################################################################################################################
#Estimate ROC(0.1)
########################################################################################################################################################################################################################
specificity_fixed=0.9

cut_off_M1=quantile(data_validation$Z1_PEB[(data_validation$D==0)],probs=specificity_fixed,type=1)
cut_off_M2=quantile(data_validation$Z2_PEB[(data_validation$D==0)],probs=specificity_fixed,type=1)
cut_off_M3=quantile(data_validation$Z3_PEB[(data_validation$D==0)],probs=specificity_fixed,type=1)

at_least_one_positive_M1=rep(NA,N-n_0)
at_least_one_positive_M2=rep(NA,N-n_0)
at_least_one_positive_M3=rep(NA,N-n_0)

for(i in 1:(N-n_0))
{
  subject_data=subset(data_validation,data_validation$ID==subj_ID[i+n_0]) 
  
  at_least_one_positive_M1[i]=as.numeric(sum(subject_data$Z1_PEB>cut_off_M1)>0)
  at_least_one_positive_M2[i]=as.numeric(sum(subject_data$Z2_PEB>cut_off_M2)>0)
  at_least_one_positive_M3[i]=as.numeric(sum(subject_data$Z3_PEB>cut_off_M3)>0)
}

round(mean(at_least_one_positive_M1)*100,2)
#73.47
round(mean(at_least_one_positive_M2)*100,2)
#63.27
round(mean(at_least_one_positive_M3)*100,2)
#63.27



