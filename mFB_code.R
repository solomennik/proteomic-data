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

set.seed(1) #set seed so each run matches reported results 

########################################################################################################################################################################################################################
#Functions for model fitting
########################################################################################################################################################################################################################
posterior_mu_theta=function(Theta,sigma_theta_2,mu_0,sigma_0_2)
{
  mu_0_star=(sigma_theta_2*mu_0)/(sigma_theta_2 + N*sigma_0_2) + (sigma_0_2*sum(Theta))/(sigma_theta_2 + N*sigma_0_2)
  sigma_0_star_2=(sigma_theta_2*sigma_0_2)/(sigma_theta_2 + N*sigma_0_2)
  list(mu_0_star=mu_0_star,sigma_0_star_2=sigma_0_star_2)
}

posterior_sigma_theta_2=function(Theta,mu_theta,a_theta,b_theta)
{
  a_theta_star=a_theta + N/2
  b_theta_star=b_theta + sum((Theta-mu_theta)^2)/2
  list(a_theta_star=a_theta_star,b_theta_star=b_theta_star)
}

posterior_sigma_2=function(Y,t,D,J,Theta,I,gamma,Tau)
{
  a_sigma=sum(J)/2
  
  temp_I=rep(c(rep(0,n_0),I),J)
  temp_gamma=rep(c(rep(0,n_0),gamma),J)
  temp_gamma[is.na(temp_gamma)]=0
  temp_Tau=rep(c(rep(0,n_0),Tau),J)
  temp_Tau[is.na(temp_Tau)]=0
  
  Theta_star=rep(Theta,J)+temp_gamma*(t-temp_Tau)*as.numeric(t>temp_Tau)
  b_sigma=sum((Y-Theta_star)^2)/2
  
  list(a_sigma=a_sigma,b_sigma=b_sigma)
}

posterior_mu_gamma=function(I,gamma,sigma_gamma_2,mu_1,sigma_1_2)
{
  n_I=sum(I)
  mu_1_star=(sigma_gamma_2*mu_1)/(sigma_gamma_2 + n_I*sigma_1_2) + (sigma_1_2*sum(log(gamma),na.rm=T))/(sigma_gamma_2 + n_I*sigma_1_2)
  sigma_1_star_2=(sigma_gamma_2*sigma_1_2)/(sigma_gamma_2 + n_I*sigma_1_2)
  list(mu_1_star=mu_1_star,sigma_1_star_2=sigma_1_star_2)
}

posterior_sigma_gamma_2=function(I,gamma,mu_gamma,a_gamma,b_gamma)
{
  n_I=sum(I)
  a_gamma_star=a_gamma+n_I/2
  b_gamma_star=b_gamma+sum((log(gamma)-mu_gamma)^2,na.rm=T)/2
  list(a_gamma_star=a_gamma_star,b_gamma_star=b_gamma_star)
}

posterior_theta=function(Y_i,t_i,D_i,J_i,mu_theta,sigma_theta_2,sigma_2,I_i,gamma_i,tau_i)
{
  Y_star_i=Y_i
  if(is.na(I_i)==FALSE & I_i==1)
  {
    Y_star_i=Y_i-gamma_i*(t_i-tau_i)*as.numeric(t_i>tau_i)
  } 
  mu_theta_star=(sigma_2*mu_theta)/(sigma_2 + J_i*sigma_theta_2) + (sigma_theta_2*sum(Y_star_i))/(sigma_2 + J_i*sigma_theta_2)
  sigma_theta_star_2=(sigma_2*sigma_theta_2)/(sigma_2 + J_i*sigma_theta_2)
  list(mu_theta_star=mu_theta_star,sigma_theta_star_2=sigma_theta_star_2)
}

log_likelihood_cases=function(Y_i,t_i,J_i,sigma_2,Theta_i,I_i,gamma_i,tau_i)
{
  log_pr=NA
  if(I_i==0)
  {
    log_pr=log((2*pi*sigma_2)^(-J_i/2)) + (-sum((Y_i-Theta_i)^2)/(2*sigma_2))
  }
  if(I_i==1)
  {
    log_pr=log((2*pi*sigma_2)^(-J_i/2)) + (-sum((Y_i-Theta_i-gamma_i*(t_i-tau_i)*as.numeric(t_i>tau_i))^2)/(2*sigma_2))
  }
  log_pr
}

update_I=function(Y_i,t_i,J_i,d_i,mu_gamma,sigma_gamma_2,mu_tau,sigma_tau_2,Tau_star,sigma_2,Theta_i,I_i,gamma_i,tau_i,Pi)
{
  #new_I=NA
  new_gamma=NA
  new_tau=NA
  if(I_i==0)
  {
    gamma_star=exp(rnorm(1,mean=mu_gamma,sd=sqrt(sigma_gamma_2)))
    tau_star=rtruncnorm(1,a=d_i-Tau_star,b=d_i,mean=d_i-mu_tau,sd=sqrt(sigma_tau_2))
    temp_log_r=log_likelihood_cases(Y_i,t_i,J_i,sigma_2,Theta_i,I_i=1,gamma_i=gamma_star,tau_i=tau_star)-log_likelihood_cases(Y_i,t_i,J_i,sigma_2,Theta_i,I_i=0,gamma_i=NA,tau_i=NA)
    log_r=min(temp_log_r+log(Pi)-log(1-Pi),log(1))
    u=runif(1,min=0,max=1)
    new_I=as.numeric(log(u)<log_r)
    if(new_I==1)
    {
      new_gamma=gamma_star
      new_tau=tau_star
    }
  }
  if(I_i==1)
  {
    temp_log_r=log_likelihood_cases(Y_i,t_i,J_i,sigma_2,Theta_i,I_i=0,gamma_i=NA,tau_i=NA)-log_likelihood_cases(Y_i,t_i,J_i,sigma_2,Theta_i,I_i=1,gamma_i=gamma_i,tau_i=tau_i)
    log_r=min(temp_log_r+log(1-Pi)-log(Pi),log(1))
    u=runif(1,min=0,max=1)
    new_I=1-as.numeric(log(u)<log_r)
    if(new_I==1)
    {
      new_gamma=gamma_i
      new_tau=tau_i
    }
  }
  list(new_I=new_I,new_gamma=new_gamma,new_tau=new_tau)
}

update_gamma=function(Y_i,t_i,J_i,mu_gamma,sigma_gamma_2,sigma_2,Theta_i,I_i,gamma_i,tau_i,delta_gamma)
{
  new_gamma=NA
  gamma_accepted=NA
  new_I=0
  new_tau=NA
  if(I_i==1)
  {
    gamma_star=exp(rnorm(1,mean=log(gamma_i),sd=delta_gamma))
    log_temp_r1=log_likelihood_cases(Y_i,t_i,J_i,sigma_2,Theta_i,I_i=1,gamma_i=gamma_star,tau_i=tau_i)-log_likelihood_cases(Y_i,t_i,J_i,sigma_2,Theta_i,I_i=1,gamma_i=gamma_i,tau_i=tau_i)
    log_temp_r2=log(dlnorm(x=gamma_star,meanlog=mu_gamma,sdlog=sqrt(sigma_gamma_2)))-log(dlnorm(x=gamma_i,meanlog=mu_gamma,sdlog=sqrt(sigma_gamma_2)))
    log_temp_r3=log(dlnorm(x=gamma_i,meanlog=log(gamma_star),sdlog=delta_gamma))-log(dlnorm(x=gamma_star,meanlog=log(gamma_i),sdlog=delta_gamma))
    log_r=min(log_temp_r1+log_temp_r2+log_temp_r3,log(1))
    u=runif(1,min=0,max=1)
    new_gamma=as.numeric(log(u)<log_r)*gamma_star + (1-as.numeric(log(u)<log_r))*gamma_i
    gamma_accepted=as.numeric(log(u)<log_r)
    new_I=I_i
    new_tau=tau_i
  }
  list(new_gamma=new_gamma,gamma_accepted=gamma_accepted,new_I=new_I,new_tau=new_tau)
}

update_tau=function(Y_i,t_i,J_i,d_i,sigma_2,mu_tau,sigma_tau_2,Theta_i,I_i,gamma_i,tau_i,delta_tau,Tau_star)
{
  new_gamma=NA
  new_I=0
  new_tau=NA
  tau_accepted=NA
  if(I_i==1)
  {
    tau_star=rtruncnorm(1,a=d_i-Tau_star,b=d_i,mean=tau_i,sd=delta_tau)
    log_temp_r1=log_likelihood_cases(Y_i,t_i,J_i,sigma_2,Theta_i,I_i=1,gamma_i=gamma_i,tau_i=tau_star)-log_likelihood_cases(Y_i,t_i,J_i,sigma_2,Theta_i,I_i=1,gamma_i=gamma_i,tau_i=tau_i)
    log_temp_r2=log(dtruncnorm(x=tau_star,a=d_i-Tau_star,b=d_i,mean=d_i-mu_tau,sd=sqrt(sigma_tau_2)))-log(dtruncnorm(x=tau_i,a=d_i-Tau_star,b=d_i,mean=d_i-mu_tau,sd=sqrt(sigma_tau_2)))
    log_temp_r3=log(dtruncnorm(x=tau_i,a=d_i-Tau_star,b=d_i,mean=tau_star,sd=delta_tau))-log(dtruncnorm(x=tau_star,a=d_i-Tau_star,b=d_i,mean=tau_i,sd=delta_tau))
    log_r=min(log_temp_r1+log_temp_r2+log_temp_r3,log(1))
    u=runif(1,min=0,max=1)
    new_tau=as.numeric(log(u)<log_r)*tau_star + (1-as.numeric(log(u)<log_r))*tau_i
    tau_accepted=as.numeric(log(u)<log_r)
    new_I=I_i
    new_gamma=gamma_i
  }
  list(new_gamma=new_gamma,tau_accepted=tau_accepted,new_I=new_I,new_tau=new_tau)
}

update_mu_tau=function(Tau,I,d,D,mu_tau,sigma_tau_2,delta_mu_tau,mu_2,sigma_2_2,Tau_star)
{
  n_I=sum(I)
  mu_tau_star=rnorm(1,mean=mu_tau,sd=delta_mu_tau)
  log_temp1=sum(log(dtruncnorm(x=Tau[I==1],a=d[D==1][I==1]-Tau_star,b=d[D==1][I==1],mean=d[D==1][I==1]-mu_tau_star,sd=sqrt(sigma_tau_2)))) + log(dnorm(x=mu_tau_star,mean=mu_2,sd=sqrt(sigma_2_2)))
  log_temp2=sum(log(dtruncnorm(x=Tau[I==1],a=d[D==1][I==1]-Tau_star,b=d[D==1][I==1],mean=d[D==1][I==1]-mu_tau,sd=sqrt(sigma_tau_2)))) + log(dnorm(x=mu_tau,mean=mu_2,sd=sqrt(sigma_2_2)))
  log_r=min(log_temp1 - log_temp2,log(1))
  u=runif(1,min=0,max=1)
  new_mu_tau=as.numeric(log(u)<log_r)*mu_tau_star + (1-as.numeric(log(u)<log_r))*mu_tau
  mu_tau_accepted=as.numeric(log(u)<log_r)
  list(new_mu_tau=new_mu_tau,mu_tau_accepted=mu_tau_accepted)
}

update_sigma_tau=function(Tau,I,d,D,mu_tau,sigma_tau_2,delta_sigma_tau,a_tau,b_tau,Tau_star)
{
  n_I=sum(I)
  sigma_tau_2_star=rtruncnorm(1,a=0,b=Inf,mean=sigma_tau_2,sd=delta_sigma_tau)
  log_temp1=sum(log(dtruncnorm(x=Tau[I==1],a=d[D==1][I==1]-Tau_star,b=d[D==1][I==1],mean=d[D==1][I==1]-mu_tau,sd=sqrt(sigma_tau_2_star)))) + log(dinvgamma(x=sigma_tau_2_star,shape=a_tau,scale=b_tau))
  log_temp2=sum(log(dtruncnorm(x=Tau[I==1],a=d[D==1][I==1]-Tau_star,b=d[D==1][I==1],mean=d[D==1][I==1]-mu_tau,sd=sqrt(sigma_tau_2)))) + log(dinvgamma(x=sigma_tau_2,shape=a_tau,scale=b_tau))
  log_temp3=log(dtruncnorm(x=sigma_tau_2,a=0,b=Inf,mean=sigma_tau_2_star,sd=delta_sigma_tau))-log(dtruncnorm(x=sigma_tau_2_star,a=0,b=Inf,mean=sigma_tau_2,sd=delta_sigma_tau))
  log_r=min(log_temp1 - log_temp2 + log_temp3,log(1))
  u=runif(1,min=0,max=1)
  new_sigma_tau_2=as.numeric(log(u)<log_r)*sigma_tau_2_star + (1-as.numeric(log(u)<log_r))*sigma_tau_2
  sigma_tau_accepted=as.numeric(log(u)<log_r)
  list(new_sigma_tau_2=new_sigma_tau_2,sigma_tau_accepted=sigma_tau_accepted)
}

#I_matrix is a matrix of indicator variables of dimension (n-n_0) X K
update_mu_I=function(I_matrix,mu_I,eta_I,delta_mu_I)
{
  mu_I_star=rnorm(1,mean=mu_I,sd=delta_mu_I)
  log_temp1=log_p_all_I_mu_eta(I_matrix,mu_I_star,eta_I) + log_p_mu_I(mu_I_star)
  log_temp2=log_p_all_I_mu_eta(I_matrix,mu_I,eta_I) + log_p_mu_I(mu_I)
  log_r=min(log_temp1 - log_temp2,log(1))
  u=runif(1,min=0,max=1)
  new_mu_I=as.numeric(log(u)<log_r)*mu_I_star + (1-as.numeric(log(u)<log_r))*mu_I
  mu_I_accepted=as.numeric(log(u)<log_r)
  list(new_mu_I=new_mu_I,mu_I_accepted=mu_I_accepted)
}

#I_matrix is a matrix of indicator variables of dimension (n-n_0) X K
log_p_all_I_mu_eta=function(I_matrix,mu_I,eta_I)
{
  temp=rep(NA,nrow(I_matrix))
  for(l in 1:nrow(I_matrix))
  {
    temp[l]=log_p_I_mu_eta(I_matrix[l,],mu_I,eta_I)
  }
  sum(temp)
}

#I_vector is a vector of indicator variables
log_p_I_mu_eta=function(I_vector,mu_I,eta_I)
{
  R=array(1,c(K,K))
  R[lower.tri(R,diag=TRUE)]=0
  
  #Normaliztion constant
  I_full=find_all_comb(K)
  c=0
  for(l in 1:nrow(I_full))
  {
    temp1=sum(I_full[l,])
    temp2=array(I_full[l,],c(K,1))
    temp3=mu_I*temp1 + eta_I*(t(temp2)%*%R%*%temp2)
    c=c+exp(temp3)    
  }
  
  #Observed data
  temp1=sum(I_vector)
  temp2=array(I_vector,c(K,1))
  mu_I*temp1 + eta_I*(t(temp2)%*%R%*%temp2) - log(c)
}

#K>1
find_all_comb=function(K)
{
  temp1=array(0,c(K+1,K))
  temp1[lower.tri(temp1,diag=FALSE)]=1
  temp2=temp1[1,]
  for(k in 2:K)
  {
    temp3=unique(permn(temp1[k,]))
    l_temp3=length(temp3)
    for(l in 1:l_temp3)
    {
      temp2=rbind(temp2,temp3[[l]])
    }
  }
  temp2=rbind(temp2,temp1[K+1,])
  array(temp2,dimnames=NULL,c(nrow(temp2),ncol(temp2)))
}

log_p_mu_I=function(mu_I)
{
  temp1=log(dbeta(x=exp(mu_I)/(1+exp(mu_I)),shape1=p_1,shape2=p_2))
  temp2=log(1+exp(mu_I))
  temp1 + mu_I - 2*temp2
}

update_eta_I=function(I_matrix,mu_I,eta_I,delta_eta_I)
{
  eta_I_star_parameters=get_beta_parameters(eta_I,delta_eta_I)
  eta_I_star=rbeta(1,shape1=eta_I_star_parameters$a,shape2=eta_I_star_parameters$b) 
  log_temp1=log_p_all_I_mu_eta(I_matrix,mu_I,eta_I_star) + log(dbeta(eta_I_star,shape1=p_3,shape2=p_4))
  log_temp2=log_p_all_I_mu_eta(I_matrix,mu_I,eta_I) + log(dbeta(eta_I,shape1=p_3,shape2=p_4))
  temp_beta_parameters1=get_beta_parameters(eta_I_star,delta_eta_I)
  temp_beta_parameters2=get_beta_parameters(eta_I,delta_eta_I)
  log_temp3=log(dbeta(eta_I,shape1=temp_beta_parameters1$a,shape2=temp_beta_parameters1$b)) - log(dbeta(eta_I_star,shape1=temp_beta_parameters2$a,shape2=temp_beta_parameters2$b))
  log_r=min(log_temp1 - log_temp2 + log_temp3,log(1))
  u=runif(1,min=0,max=1)
  new_eta_I=as.numeric(log(u)<log_r)*eta_I_star + (1-as.numeric(log(u)<log_r))*eta_I
  eta_I_accepted=as.numeric(log(u)<log_r)
  list(new_eta_I=new_eta_I,eta_I_accepted=eta_I_accepted)
}

get_beta_parameters=function(mean_eta_I,std_dev_eta_I)
{
  var_eta_I=std_dev_eta_I^2
  a=mean_eta_I + (mean_eta_I^2)/var_eta_I - (mean_eta_I^3)/var_eta_I
  b=a*(1/mean_eta_I - 1)
  list(a=a, b=b)
}

get_p1=function(mu_I,eta_I)
{
  c_inv=1+3*exp(mu_I)+3*exp(2*mu_I+eta_I)+exp(3*mu_I+3*eta_I)
  (exp(mu_I)+2*exp(2*mu_I+eta_I)+exp(3*mu_I+3*eta_I))/c_inv
}

get_p2_I_1=function(I_1,mu_I,eta_I)
{
  c_inv=1+3*exp(mu_I)+3*exp(2*mu_I+eta_I)+exp(3*mu_I+3*eta_I)
  p_I_1=(exp(mu_I)+2*exp(2*mu_I+eta_I)+exp(3*mu_I+3*eta_I))/c_inv
  num=(exp(mu_I*(I_1+1)+eta_I*I_1) + exp(mu_I*(I_1+2)+eta_I*(2*I_1 + 1)))/c_inv
  dem=p_I_1*I_1 + (1-p_I_1)*(1-I_1)
  num/dem
}

get_p3_I_1_2=function(I_1,I_2,mu_I,eta_I)
{
  F_I_3=mu_I + eta_I*(I_1 + I_2)
  exp(F_I_3)/(1+exp(F_I_3))
}

########################################################################################################################################################################################################################
#Analyze Training Data (pg 6 Supplementary Materials)
########################################################################################################################################################################################################################
#Get sample size 
N=length(unique(data$ID)) #total number of patients
n_0=length(unique(data$ID[(data$D==0)])) #number of control patients

#Extract vectors
D=data$D[(data$obs_number==1)] #case/control status
d=data$d[(data$obs_number==1)] #diagnosis/follow-up times
J=data$J[(data$obs_number==1)] #number of screening visits
K=3 #Number of biomarkers in the study

#Specify fixed parameters
p_1=30
p_2=30
p_3=5
p_4=45

#Biomarker 1
mu_0_1=2
sigma_0_1_2=0.1
a_theta_1=2
b_theta_1=0.05
mu_1_1=log(16)
sigma_1_1_2=0.1 
a_gamma_1=4
b_gamma_1=0.1
Tau_star_1=2
mu_2_1=1
sigma_2_1_2=0.1
a_tau_1=10
b_tau_1=(0.75^2)*(a_tau_1+1)

#Biomarker 2
mu_0_2=3
sigma_0_2_2=0.1
a_theta_2=2
b_theta_2=0.05
mu_1_2=log(16)
sigma_1_2_2=0.05
a_gamma_2=2*2
b_gamma_2=0.05*2
Tau_star_2=2
mu_2_2=1
sigma_2_2_2=0.1
a_tau_2=10
b_tau_2=(0.75^2)*(a_tau_2+1)

#Biomarker 3
mu_0_3=2
sigma_0_3_2=0.1
a_theta_3=2
b_theta_3=0.05
mu_1_3=log(16)
sigma_1_3_2=0.05
a_gamma_3=2*2
b_gamma_3=0.05*2
Tau_star_3=2
mu_2_3=1
sigma_2_3_2=0.1
a_tau_3=10
b_tau_3=(0.75^2)*(a_tau_1+1)

S=10000
delta_gamma_set=0.1
delta_tau_set=0.5
delta_mu_tau_set=0.75
delta_sigma_tau_set=0.3
delta_mu_I_set=0.2
delta_eta_I_set=0.02

###########################Initialize parameters from priors for first chain
temp_p=rbeta(n=1,shape1=p_1,shape2=p_2)
mu_I_out=log(temp_p/(1-temp_p))
eta_I_out=rbeta(n=1,shape1=p_3,shape2=p_4)

#biomarker 1
mu_theta_1_out=rnorm(1,mean=mu_0_1,sd=sqrt(sigma_0_1_2))
sigma_theta_1_2_out=rinvgamma(1,shape=a_theta_1,scale=b_theta_1)
Theta_1_out=array(rnorm(N,mean=mu_theta_1_out,sd=sigma_theta_1_2_out),c(N,1))
I_1_out=array(rbinom(3*(N-n_0),size=1,prob=0.5),c(N-n_0,3))
mu_gamma_1_out=rnorm(1,mean=mu_1_1,sd=sqrt(sigma_1_1_2))
sigma_gamma_1_2_out=rinvgamma(1,shape=a_gamma_1,scale=b_gamma_1)
mu_tau_1_out=1
sigma_tau_1_2_out=(0.75)^2
gamma_1_out=array(exp(rnorm(3*(N-n_0),mean=mu_gamma_1_out,sd=sqrt(sigma_gamma_1_2_out))),c(N-n_0,3))
gamma_1_out[I_1_out[,1]==0,1]=NA
gamma_1_out[I_1_out[,2]==0,2]=NA
gamma_1_out[I_1_out[,3]==0,3]=NA
Tau_1_out=array(rtruncnorm(n=3*(N-n_0),a=d[D==1]-Tau_star_1,b=d[D==1],mean=d[D==1]-mu_tau_1_out,sd=sqrt(sigma_tau_1_2_out)),c(N-n_0,3))
Tau_1_out[I_1_out[,1]==0,1]=NA
Tau_1_out[I_1_out[,2]==0,2]=NA
Tau_1_out[I_1_out[,3]==0,3]=NA
sigma_1_2_out=0.1

#biomarker 2
mu_theta_2_out=rnorm(1,mean=mu_0_2,sd=sqrt(sigma_0_2_2))
sigma_theta_2_2_out=rinvgamma(1,shape=a_theta_2,scale=b_theta_2)
Theta_2_out=array(rnorm(N,mean=mu_theta_2_out,sd=sigma_theta_2_2_out),c(N,1))
I_2_out=array(rbinom(3*(N-n_0),size=1,prob=0.5),c(N-n_0,3))
mu_gamma_2_out=rnorm(1,mean=mu_1_2,sd=sqrt(sigma_1_2_2))
sigma_gamma_2_2_out=rinvgamma(1,shape=a_gamma_2,scale=b_gamma_2)
mu_tau_2_out=1
sigma_tau_2_2_out=(0.75)^2
gamma_2_out=array(exp(rnorm(3*(N-n_0),mean=mu_gamma_2_out,sd=sqrt(sigma_gamma_2_2_out))),c(N-n_0,3))
gamma_2_out[I_2_out[,1]==0,1]=NA
gamma_2_out[I_2_out[,2]==0,2]=NA
gamma_2_out[I_2_out[,3]==0,3]=NA
Tau_2_out=array(rtruncnorm(n=3*(N-n_0),a=d[D==1]-Tau_star_1,b=d[D==1],mean=d[D==1]-mu_tau_2_out,sd=sqrt(sigma_tau_2_2_out)),c(N-n_0,3))
Tau_2_out[I_2_out[,1]==0,1]=NA
Tau_2_out[I_2_out[,2]==0,2]=NA
Tau_2_out[I_2_out[,3]==0,3]=NA
sigma_2_2_out=0.1

#biomarker 3
mu_theta_3_out=rnorm(1,mean=mu_0_3,sd=sqrt(sigma_0_3_2))
sigma_theta_3_2_out=rinvgamma(1,shape=a_theta_3,scale=b_theta_3)
Theta_3_out=array(rnorm(N,mean=mu_theta_3_out,sd=sigma_theta_3_2_out),c(N,1))
I_3_out=array(rbinom(3*(N-n_0),size=1,prob=0.5),c(N-n_0,3))
mu_gamma_3_out=rnorm(1,mean=mu_1_3,sd=sqrt(sigma_1_3_2))
sigma_gamma_3_2_out=rinvgamma(1,shape=a_gamma_3,scale=b_gamma_3)
mu_tau_3_out=1
sigma_tau_3_2_out=(0.75)^2
gamma_3_out=array(exp(rnorm(3*(N-n_0),mean=mu_gamma_3_out,sd=sqrt(sigma_gamma_3_2_out))),c(N-n_0,3))
gamma_3_out[I_3_out[,1]==0,1]=NA
gamma_3_out[I_3_out[,2]==0,2]=NA
gamma_3_out[I_3_out[,3]==0,3]=NA
Tau_3_out=array(rtruncnorm(n=3*(N-n_0),a=d[D==1]-Tau_star_3,b=d[D==1],mean=d[D==1]-mu_tau_3_out,sd=sqrt(sigma_tau_3_2_out)),c(N-n_0,3))
Tau_3_out[I_3_out[,1]==0,1]=NA
Tau_3_out[I_3_out[,2]==0,2]=NA
Tau_3_out[I_3_out[,3]==0,3]=NA
sigma_3_2_out=0.1

#Update parameters
for(i in 1:S)
{
  #MU_THETA_1
  mu_theta_parameters=posterior_mu_theta(Theta=Theta_1_out[,i],sigma_theta_2=sigma_theta_1_2_out[i],mu_0=mu_0_1,sigma_0_2=sigma_0_1_2)
  mu_theta_1_out=c(mu_theta_1_out,rnorm(1,mean=mu_theta_parameters$mu_0_star,sd=sqrt(mu_theta_parameters$sigma_0_star_2)))
  
  #MU_THETA_2
  mu_theta_parameters=posterior_mu_theta(Theta=Theta_2_out[,i],sigma_theta_2=sigma_theta_2_2_out[i],mu_0=mu_0_2,sigma_0_2=sigma_0_2_2)
  mu_theta_2_out=c(mu_theta_2_out,rnorm(1,mean=mu_theta_parameters$mu_0_star,sd=sqrt(mu_theta_parameters$sigma_0_star_2)))
  
  #MU_THETA_3
  mu_theta_parameters=posterior_mu_theta(Theta=Theta_3_out[,i],sigma_theta_2=sigma_theta_3_2_out[i],mu_0=mu_0_3,sigma_0_2=sigma_0_3_2)
  mu_theta_3_out=c(mu_theta_3_out,rnorm(1,mean=mu_theta_parameters$mu_0_star,sd=sqrt(mu_theta_parameters$sigma_0_star_2)))
  
  #SIGMA_THETA_1
  sigma_theta_parameters=posterior_sigma_theta_2(Theta=Theta_1_out[,i],mu_theta=mu_theta_1_out[i+1],a_theta=a_theta_1,b_theta=b_theta_1)
  sigma_theta_1_2_out=c(sigma_theta_1_2_out,rinvgamma(1,shape=sigma_theta_parameters$a_theta_star,scale=sigma_theta_parameters$b_theta_star))
  
  #SIGMA_THETA_2
  sigma_theta_parameters=posterior_sigma_theta_2(Theta=Theta_2_out[,i],mu_theta=mu_theta_2_out[i+1],a_theta=a_theta_2,b_theta=b_theta_2)
  sigma_theta_2_2_out=c(sigma_theta_2_2_out,rinvgamma(1,shape=sigma_theta_parameters$a_theta_star,scale=sigma_theta_parameters$b_theta_star))
  
  #SIGMA_THETA_3
  sigma_theta_parameters=posterior_sigma_theta_2(Theta=Theta_3_out[,i],mu_theta=mu_theta_3_out[i+1],a_theta=a_theta_3,b_theta=b_theta_3)
  sigma_theta_3_2_out=c(sigma_theta_3_2_out,rinvgamma(1,shape=sigma_theta_parameters$a_theta_star,scale=sigma_theta_parameters$b_theta_star))
  
  #SIGMA_1
  sigma_parameters=posterior_sigma_2(Y=data$Y1,t=data$t,D=D,J=J,Theta=Theta_1_out[,i],I=I_1_out[,3*i],gamma=gamma_1_out[,3*i],Tau=Tau_1_out[,3*i])
  sigma_1_2_out=c(sigma_1_2_out,rinvgamma(1,shape=sigma_parameters$a_sigma,scale=sigma_parameters$b_sigma))
  
  #SIGMA_2
  sigma_parameters=posterior_sigma_2(Y=data$Y2,t=data$t,D=D,J=J,Theta=Theta_2_out[,i],I=I_2_out[,3*i],gamma=gamma_2_out[,3*i],Tau=Tau_2_out[,3*i])
  sigma_2_2_out=c(sigma_2_2_out,rinvgamma(1,shape=sigma_parameters$a_sigma,scale=sigma_parameters$b_sigma))
  
  #SIGMA_3
  sigma_parameters=posterior_sigma_2(Y=data$Y3,t=data$t,D=D,J=J,Theta=Theta_3_out[,i],I=I_3_out[,3*i],gamma=gamma_3_out[,3*i],Tau=Tau_3_out[,3*i])
  sigma_3_2_out=c(sigma_3_2_out,rinvgamma(1,shape=sigma_parameters$a_sigma,scale=sigma_parameters$b_sigma))
  
  #MU_I
  temp_I_matrix=cbind(I_1_out[,3*i],I_2_out[,3*i],I_3_out[,3*i])
  mu_I_output=update_mu_I(temp_I_matrix,mu_I_out[i],eta_I_out[i],delta_mu_I_set)
  mu_I_out=c(mu_I_out,mu_I_output$new_mu_I)
  
  #ETA_I
  temp_I_matrix=cbind(I_1_out[,3*i],I_2_out[,3*i],I_3_out[,3*i])
  eta_I_output=update_eta_I(temp_I_matrix,mu_I_out[i+1],eta_I_out[i],delta_eta_I_set)
  eta_I_out=c(eta_I_out,eta_I_output$new_eta_I)
  
  #MU_GAMMA_1
  mu_gamma_parameters=posterior_mu_gamma(I=I_1_out[,3*i],gamma=gamma_1_out[,3*i],sigma_gamma_2=sigma_gamma_1_2_out[i],mu_1=mu_1_1,sigma_1_2=sigma_1_1_2)
  mu_gamma_1_out=c(mu_gamma_1_out,rnorm(1,mean=mu_gamma_parameters$mu_1_star,sd=sqrt(mu_gamma_parameters$sigma_1_star_2)))
  
  #MU_GAMMA_2
  mu_gamma_parameters=posterior_mu_gamma(I=I_2_out[,3*i],gamma=gamma_2_out[,3*i],sigma_gamma_2=sigma_gamma_2_2_out[i],mu_1=mu_1_2,sigma_1_2=sigma_1_2_2)
  mu_gamma_2_out=c(mu_gamma_2_out,rnorm(1,mean=mu_gamma_parameters$mu_1_star,sd=sqrt(mu_gamma_parameters$sigma_1_star_2)))
  
  #MU_GAMMA_3
  mu_gamma_parameters=posterior_mu_gamma(I=I_3_out[,3*i],gamma=gamma_3_out[,3*i],sigma_gamma_2=sigma_gamma_3_2_out[i],mu_1=mu_1_3,sigma_1_2=sigma_1_3_2)
  mu_gamma_3_out=c(mu_gamma_3_out,rnorm(1,mean=mu_gamma_parameters$mu_1_star,sd=sqrt(mu_gamma_parameters$sigma_1_star_2)))
  
  #SIGMA_GAMMA_1
  sigma_gamma_parameters=posterior_sigma_gamma_2(I=I_1_out[,3*i],gamma=gamma_1_out[,3*i],mu_gamma=mu_gamma_1_out[i+1],a_gamma=a_gamma_1,b_gamma=b_gamma_1)
  sigma_gamma_1_2_out=c(sigma_gamma_1_2_out,rinvgamma(1,shape=sigma_gamma_parameters$a_gamma_star,scale=sigma_gamma_parameters$b_gamma_star))
  
  #SIGMA_GAMMA_2
  sigma_gamma_parameters=posterior_sigma_gamma_2(I=I_2_out[,3*i],gamma=gamma_2_out[,3*i],mu_gamma=mu_gamma_2_out[i+1],a_gamma=a_gamma_2,b_gamma=b_gamma_2)
  sigma_gamma_2_2_out=c(sigma_gamma_2_2_out,rinvgamma(1,shape=sigma_gamma_parameters$a_gamma_star,scale=sigma_gamma_parameters$b_gamma_star))
  
  #SIGMA_GAMMA_3
  sigma_gamma_parameters=posterior_sigma_gamma_2(I=I_3_out[,3*i],gamma=gamma_3_out[,3*i],mu_gamma=mu_gamma_3_out[i+1],a_gamma=a_gamma_3,b_gamma=b_gamma_3)
  sigma_gamma_3_2_out=c(sigma_gamma_3_2_out,rinvgamma(1,shape=sigma_gamma_parameters$a_gamma_star,scale=sigma_gamma_parameters$b_gamma_star))
  
  #MU_TAU_1
  mu_tau_output=update_mu_tau(Tau=Tau_1_out[,3*i],I=I_1_out[,3*i],d=d,D=D,mu_tau=mu_tau_1_out[i],sigma_tau_2=sigma_tau_1_2_out[i],delta_mu_tau=delta_mu_tau_set,mu_2=mu_2_1,sigma_2_2=sigma_2_1_2,Tau_star=Tau_star_1)
  mu_tau_1_out=c(mu_tau_1_out,mu_tau_output$new_mu_tau)
  
  #MU_TAU_2
  mu_tau_output=update_mu_tau(Tau=Tau_2_out[,3*i],I=I_2_out[,3*i],d=d,D=D,mu_tau=mu_tau_2_out[i],sigma_tau_2=sigma_tau_2_2_out[i],delta_mu_tau=delta_mu_tau_set,mu_2=mu_2_2,sigma_2_2=sigma_2_2_2,Tau_star=Tau_star_2)
  mu_tau_2_out=c(mu_tau_2_out,mu_tau_output$new_mu_tau)
  
  #MU_TAU_3
  mu_tau_output=update_mu_tau(Tau=Tau_3_out[,3*i],I=I_3_out[,3*i],d=d,D=D,mu_tau=mu_tau_3_out[i],sigma_tau_2=sigma_tau_3_2_out[i],delta_mu_tau=delta_mu_tau_set,mu_2=mu_2_3,sigma_2_2=sigma_2_3_2,Tau_star=Tau_star_3)
  mu_tau_3_out=c(mu_tau_3_out,mu_tau_output$new_mu_tau)
  
  #SIGMA_TAU_1
  sigma_tau_output=update_sigma_tau(Tau=Tau_1_out[,3*i],I=I_1_out[,3*i],d=d,D=D,mu_tau=mu_tau_1_out[i+1],sigma_tau_2=sigma_tau_1_2_out[i],delta_sigma_tau=delta_sigma_tau_set,a_tau=a_tau_1,b_tau=b_tau_1,Tau_star=Tau_star_1)
  sigma_tau_1_2_out=c(sigma_tau_1_2_out,sigma_tau_output$new_sigma_tau_2)
  
  #SIGMA_TAU_2
  sigma_tau_output=update_sigma_tau(Tau=Tau_2_out[,3*i],I=I_2_out[,3*i],d=d,D=D,mu_tau=mu_tau_2_out[i+1],sigma_tau_2=sigma_tau_2_2_out[i],delta_sigma_tau=delta_sigma_tau_set,a_tau=a_tau_2,b_tau=b_tau_2,Tau_star=Tau_star_2)
  sigma_tau_2_2_out=c(sigma_tau_2_2_out,sigma_tau_output$new_sigma_tau_2)
  
  #SIGMA_TAU_3
  sigma_tau_output=update_sigma_tau(Tau=Tau_3_out[,3*i],I=I_3_out[,3*i],d=d,D=D,mu_tau=mu_tau_3_out[i+1],sigma_tau_2=sigma_tau_3_2_out[i],delta_sigma_tau=delta_sigma_tau_set,a_tau=a_tau_3,b_tau=b_tau_3,Tau_star=Tau_star_3)
  sigma_tau_3_2_out=c(sigma_tau_3_2_out,sigma_tau_output$new_sigma_tau_2)
  
  Theta_1_new=rep(NA,N)
  Theta_2_new=rep(NA,N)
  Theta_3_new=rep(NA,N)
  
  I_1_new_1=rep(NA,N-n_0)
  I_1_new_2=rep(NA,N-n_0)
  I_1_new_3=rep(NA,N-n_0)
  
  I_2_new_1=rep(NA,N-n_0)
  I_2_new_2=rep(NA,N-n_0)
  I_2_new_3=rep(NA,N-n_0)
  
  I_3_new_1=rep(NA,N-n_0)
  I_3_new_2=rep(NA,N-n_0)
  I_3_new_3=rep(NA,N-n_0)
  
  gamma_1_new_1=rep(NA,N-n_0)
  gamma_1_new_2=rep(NA,N-n_0)
  gamma_1_new_3=rep(NA,N-n_0)
  
  gamma_2_new_1=rep(NA,N-n_0)
  gamma_2_new_2=rep(NA,N-n_0)
  gamma_2_new_3=rep(NA,N-n_0)
  
  gamma_3_new_1=rep(NA,N-n_0)
  gamma_3_new_2=rep(NA,N-n_0)
  gamma_3_new_3=rep(NA,N-n_0)
  
  Tau_1_new_1=rep(NA,N-n_0)
  Tau_1_new_2=rep(NA,N-n_0)
  Tau_1_new_3=rep(NA,N-n_0)
  
  Tau_2_new_1=rep(NA,N-n_0)
  Tau_2_new_2=rep(NA,N-n_0)
  Tau_2_new_3=rep(NA,N-n_0)
  
  Tau_3_new_1=rep(NA,N-n_0)
  Tau_3_new_2=rep(NA,N-n_0)
  Tau_3_new_3=rep(NA,N-n_0)
  
  for(j in 1:n_0)
  {
    #THETA_1_j
    theta_j_parameters=posterior_theta(Y_i=data$Y1[(data$ID==j)],t_i=data$t[(data$ID==j)],D_i=D[j],J_i=J[j],mu_theta=mu_theta_1_out[i+1],sigma_theta_2=sigma_theta_1_2_out[i+1],sigma_2=sigma_1_2_out[i+1],I_i=NA,gamma_i=NA,tau_i=NA)
    Theta_1_new[j]=rnorm(1,mean=theta_j_parameters$mu_theta_star,sd=sqrt(theta_j_parameters$sigma_theta_star_2))
    
    #THETA_2_j
    theta_j_parameters=posterior_theta(Y_i=data$Y2[(data$ID==j)],t_i=data$t[(data$ID==j)],D_i=D[j],J_i=J[j],mu_theta=mu_theta_2_out[i+1],sigma_theta_2=sigma_theta_2_2_out[i+1],sigma_2=sigma_2_2_out[i+1],I_i=NA,gamma_i=NA,tau_i=NA)
    Theta_2_new[j]=rnorm(1,mean=theta_j_parameters$mu_theta_star,sd=sqrt(theta_j_parameters$sigma_theta_star_2))
    
    #THETA_3_j
    theta_j_parameters=posterior_theta(Y_i=data$Y3[(data$ID==j)],t_i=data$t[(data$ID==j)],D_i=D[j],J_i=J[j],mu_theta=mu_theta_3_out[i+1],sigma_theta_2=sigma_theta_3_2_out[i+1],sigma_2=sigma_3_2_out[i+1],I_i=NA,gamma_i=NA,tau_i=NA)
    Theta_3_new[j]=rnorm(1,mean=theta_j_parameters$mu_theta_star,sd=sqrt(theta_j_parameters$sigma_theta_star_2))
  }
  for(j in (n_0+1):N)
  {
    #THETA_1_j
    theta_j_parameters=posterior_theta(Y_i=data$Y1[(data$ID==j)],t_i=data$t[(data$ID==j)],D_i=D[j],J_i=J[j],mu_theta=mu_theta_1_out[i+1],sigma_theta_2=sigma_theta_1_2_out[i+1],sigma_2=sigma_1_2_out[i+1],I_i=I_1_out[j-n_0,3*i],gamma_i=gamma_1_out[j-n_0,3*i],tau_i=Tau_1_out[j-n_0,3*i])
    Theta_1_new[j]=rnorm(1,mean=theta_j_parameters$mu_theta_star,sd=sqrt(theta_j_parameters$sigma_theta_star_2))
    
    #THETA_2_j
    theta_j_parameters=posterior_theta(Y_i=data$Y2[(data$ID==j)],t_i=data$t[(data$ID==j)],D_i=D[j],J_i=J[j],mu_theta=mu_theta_2_out[i+1],sigma_theta_2=sigma_theta_2_2_out[i+1],sigma_2=sigma_2_2_out[i+1],I_i=I_2_out[j-n_0,3*i],gamma_i=gamma_2_out[j-n_0,3*i],tau_i=Tau_2_out[j-n_0,3*i])
    Theta_2_new[j]=rnorm(1,mean=theta_j_parameters$mu_theta_star,sd=sqrt(theta_j_parameters$sigma_theta_star_2))
    
    #THETA_3_j
    theta_j_parameters=posterior_theta(Y_i=data$Y3[(data$ID==j)],t_i=data$t[(data$ID==j)],D_i=D[j],J_i=J[j],mu_theta=mu_theta_3_out[i+1],sigma_theta_2=sigma_theta_3_2_out[i+1],sigma_2=sigma_3_2_out[i+1],I_i=I_3_out[j-n_0,3*i],gamma_i=gamma_3_out[j-n_0,3*i],tau_i=Tau_3_out[j-n_0,3*i])
    Theta_3_new[j]=rnorm(1,mean=theta_j_parameters$mu_theta_star,sd=sqrt(theta_j_parameters$sigma_theta_star_2))
    
    #I_1
    F_temp=mu_I_out[i+1] + eta_I_out[i+1]*(I_2_out[j-n_0,3*i] + I_3_out[j-n_0,3*i])
    Pi_1_input=exp(F_temp)/(1+exp(F_temp))
    I_output=update_I(Y_i=data$Y1[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],d_i=d[j],mu_gamma=mu_gamma_1_out[i+1],sigma_gamma_2=sigma_gamma_1_2_out[i+1],mu_tau=mu_tau_1_out[i+1],sigma_tau_2=sigma_tau_1_2_out[i+1],Tau_star=Tau_star_1,sigma_2=sigma_1_2_out[i+1],Theta_i=Theta_1_new[j],I_i=I_1_out[j-n_0,3*i],gamma_i=gamma_1_out[j-n_0,3*i],tau_i=Tau_1_out[j-n_0,3*i],Pi=Pi_1_input)
    I_1_new_1[j-n_0]=I_output$new_I
    Tau_1_new_1[j-n_0]=I_output$new_tau
    gamma_1_new_1[j-n_0]=I_output$new_gamma
    
    #I_2
    F_temp=mu_I_out[i+1] + eta_I_out[i+1]*(I_1_new_1[j-n_0] + I_3_out[j-n_0,3*i])
    Pi_2_input=exp(F_temp)/(1+exp(F_temp))
    I_output=update_I(Y_i=data$Y2[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],d_i=d[j],mu_gamma=mu_gamma_2_out[i+1],sigma_gamma_2=sigma_gamma_2_2_out[i+1],mu_tau=mu_tau_2_out[i+1],sigma_tau_2=sigma_tau_2_2_out[i+1],Tau_star=Tau_star_2,sigma_2=sigma_2_2_out[i+1],Theta_i=Theta_2_new[j],I_i=I_2_out[j-n_0,3*i],gamma_i=gamma_2_out[j-n_0,3*i],tau_i=Tau_2_out[j-n_0,3*i],Pi=Pi_2_input)
    I_2_new_1[j-n_0]=I_output$new_I
    Tau_2_new_1[j-n_0]=I_output$new_tau
    gamma_2_new_1[j-n_0]=I_output$new_gamma
    
    #I_3
    F_temp=mu_I_out[i+1] + eta_I_out[i+1]*(I_1_new_1[j-n_0] + I_2_new_1[j-n_0])
    Pi_3_input=exp(F_temp)/(1+exp(F_temp))
    I_output=update_I(Y_i=data$Y3[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],d_i=d[j],mu_gamma=mu_gamma_3_out[i+1],sigma_gamma_2=sigma_gamma_3_2_out[i+1],mu_tau=mu_tau_3_out[i+1],sigma_tau_2=sigma_tau_3_2_out[i+1],Tau_star=Tau_star_3,sigma_2=sigma_3_2_out[i+1],Theta_i=Theta_3_new[j],I_i=I_3_out[j-n_0,3*i],gamma_i=gamma_3_out[j-n_0,3*i],tau_i=Tau_3_out[j-n_0,3*i],Pi=Pi_3_input)
    I_3_new_1[j-n_0]=I_output$new_I
    Tau_3_new_1[j-n_0]=I_output$new_tau
    gamma_3_new_1[j-n_0]=I_output$new_gamma
    
    #GAMMA_1
    gamma_output=update_gamma(Y_i=data$Y1[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],mu_gamma=mu_gamma_1_out[i+1],sigma_gamma_2=sigma_gamma_1_2_out[i+1],sigma_2=sigma_1_2_out[i+1],Theta_i=Theta_1_new[j],I_i=I_1_new_1[j-n_0],gamma_i=gamma_1_new_1[j-n_0],tau_i=Tau_1_new_1[j-n_0],delta_gamma=delta_gamma_set)
    gamma_1_new_2[j-n_0]=gamma_output$new_gamma
    I_1_new_2[j-n_0]=gamma_output$new_I
    Tau_1_new_2[j-n_0]=gamma_output$new_tau
    
    #GAMMA_2
    gamma_output=update_gamma(Y_i=data$Y2[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],mu_gamma=mu_gamma_2_out[i+1],sigma_gamma_2=sigma_gamma_2_2_out[i+1],sigma_2=sigma_2_2_out[i+1],Theta_i=Theta_2_new[j],I_i=I_2_new_1[j-n_0],gamma_i=gamma_2_new_1[j-n_0],tau_i=Tau_2_new_1[j-n_0],delta_gamma=delta_gamma_set)
    gamma_2_new_2[j-n_0]=gamma_output$new_gamma
    I_2_new_2[j-n_0]=gamma_output$new_I
    Tau_2_new_2[j-n_0]=gamma_output$new_tau
    
    #GAMMA_3
    gamma_output=update_gamma(Y_i=data$Y3[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],mu_gamma=mu_gamma_3_out[i+1],sigma_gamma_2=sigma_gamma_3_2_out[i+1],sigma_2=sigma_3_2_out[i+1],Theta_i=Theta_3_new[j],I_i=I_3_new_1[j-n_0],gamma_i=gamma_3_new_1[j-n_0],tau_i=Tau_3_new_1[j-n_0],delta_gamma=delta_gamma_set)
    gamma_3_new_2[j-n_0]=gamma_output$new_gamma
    I_3_new_2[j-n_0]=gamma_output$new_I
    Tau_3_new_2[j-n_0]=gamma_output$new_tau
    
    #TAU_1
    tau_output=update_tau(Y_i=data$Y1[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],d_i=d[j],sigma_2=sigma_1_2_out[i+1],mu_tau=mu_tau_1_out[i+1],sigma_tau_2=sigma_tau_1_2_out[i+1],Theta_i=Theta_1_new[j],I_i=I_1_new_2[j-n_0],gamma_i=gamma_1_new_2[j-n_0],tau_i=Tau_1_new_2[j-n_0],delta_tau=delta_tau_set,Tau_star=Tau_star_1)
    I_1_new_3[j-n_0]=tau_output$new_I
    Tau_1_new_3[j-n_0]=tau_output$new_tau
    gamma_1_new_3[j-n_0]=tau_output$new_gamma
    
    #TAU_2
    tau_output=update_tau(Y_i=data$Y2[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],d_i=d[j],sigma_2=sigma_2_2_out[i+1],mu_tau=mu_tau_2_out[i+1],sigma_tau_2=sigma_tau_2_2_out[i+1],Theta_i=Theta_2_new[j],I_i=I_2_new_2[j-n_0],gamma_i=gamma_2_new_2[j-n_0],tau_i=Tau_2_new_2[j-n_0],delta_tau=delta_tau_set,Tau_star=Tau_star_2)
    I_2_new_3[j-n_0]=tau_output$new_I
    Tau_2_new_3[j-n_0]=tau_output$new_tau
    gamma_2_new_3[j-n_0]=tau_output$new_gamma
    
    #TAU_1
    tau_output=update_tau(Y_i=data$Y3[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],d_i=d[j],sigma_2=sigma_3_2_out[i+1],mu_tau=mu_tau_3_out[i+1],sigma_tau_2=sigma_tau_3_2_out[i+1],Theta_i=Theta_3_new[j],I_i=I_3_new_2[j-n_0],gamma_i=gamma_3_new_2[j-n_0],tau_i=Tau_3_new_2[j-n_0],delta_tau=delta_tau_set,Tau_star=Tau_star_3)
    I_3_new_3[j-n_0]=tau_output$new_I
    Tau_3_new_3[j-n_0]=tau_output$new_tau
    gamma_3_new_3[j-n_0]=tau_output$new_gamma
    
  }
  Theta_1_out=cbind(Theta_1_out,Theta_1_new)
  Theta_2_out=cbind(Theta_2_out,Theta_2_new)
  Theta_3_out=cbind(Theta_3_out,Theta_3_new)
  
  I_1_out=cbind(I_1_out,I_1_new_1,I_1_new_2,I_1_new_3)
  I_2_out=cbind(I_2_out,I_2_new_1,I_2_new_2,I_2_new_3)
  I_3_out=cbind(I_3_out,I_3_new_1,I_3_new_2,I_3_new_3)
  
  gamma_1_out=cbind(gamma_1_out,gamma_1_new_1,gamma_1_new_2,gamma_1_new_3)
  gamma_2_out=cbind(gamma_2_out,gamma_2_new_1,gamma_2_new_2,gamma_2_new_3)
  gamma_3_out=cbind(gamma_3_out,gamma_3_new_1,gamma_3_new_2,gamma_3_new_3)
  
  Tau_1_out=cbind(Tau_1_out,Tau_1_new_1,Tau_1_new_2,Tau_1_new_3)
  Tau_2_out=cbind(Tau_2_out,Tau_2_new_1,Tau_2_new_2,Tau_2_new_3)
  Tau_3_out=cbind(Tau_3_out,Tau_3_new_1,Tau_3_new_2,Tau_3_new_3)
  
  print(i)
}

###########################Initialize parameters from priors for second chain
temp_p=rbeta(n=1,shape1=p_1,shape2=p_2)
mu_I_out_2=log(temp_p/(1-temp_p))
eta_I_out_2=rbeta(n=1,shape1=p_3,shape2=p_4)

#Biomarker 1
mu_theta_1_out_2=rnorm(1,mean=mu_0_1,sd=sqrt(sigma_0_1_2))
sigma_theta_1_2_out_2=rinvgamma(1,shape=a_theta_1,scale=b_theta_1)
Theta_1_out_2=array(rnorm(N,mean=mu_theta_1_out_2,sd=sigma_theta_1_2_out_2),c(N,1))
I_1_out_2=array(rbinom(3*(N-n_0),size=1,prob=0.5),c(N-n_0,3))
mu_gamma_1_out_2=rnorm(1,mean=mu_1_1,sd=sqrt(sigma_1_1_2))
sigma_gamma_1_2_out_2=rinvgamma(1,shape=a_gamma_1,scale=b_gamma_1)
mu_tau_1_out_2=2
sigma_tau_1_2_out_2=(2)^2
gamma_1_out_2=array(exp(rnorm(3*(N-n_0),mean=mu_gamma_1_out_2,sd=sqrt(sigma_gamma_1_2_out_2))),c(N-n_0,3))
gamma_1_out_2[I_1_out_2[,1]==0,1]=NA
gamma_1_out_2[I_1_out_2[,2]==0,2]=NA
gamma_1_out_2[I_1_out_2[,3]==0,3]=NA
Tau_1_out_2=array(rtruncnorm(n=3*(N-n_0),a=d[D==1]-Tau_star_1,b=d[D==1],mean=d[D==1]-mu_tau_1_out_2,sd=sqrt(sigma_tau_1_2_out_2)),c(N-n_0,3))
Tau_1_out_2[I_1_out_2[,1]==0,1]=NA
Tau_1_out_2[I_1_out_2[,2]==0,2]=NA
Tau_1_out_2[I_1_out_2[,3]==0,3]=NA
sigma_1_2_out_2=1

#Biomarker 2
mu_theta_2_out_2=rnorm(1,mean=mu_0_2,sd=sqrt(sigma_0_2_2))
sigma_theta_2_2_out_2=rinvgamma(1,shape=a_theta_2,scale=b_theta_2)
Theta_2_out_2=array(rnorm(N,mean=mu_theta_2_out_2,sd=sigma_theta_2_2_out_2),c(N,1))
I_2_out_2=array(rbinom(3*(N-n_0),size=1,prob=0.5),c(N-n_0,3))
mu_gamma_2_out_2=rnorm(1,mean=mu_1_2,sd=sqrt(sigma_1_2_2))
sigma_gamma_2_2_out_2=rinvgamma(1,shape=a_gamma_2,scale=b_gamma_2)
mu_tau_2_out_2=0.5
sigma_tau_2_2_out_2=1
gamma_2_out_2=array(exp(rnorm(3*(N-n_0),mean=mu_gamma_2_out_2,sd=sqrt(sigma_gamma_2_2_out_2))),c(N-n_0,3))
gamma_2_out_2[I_2_out_2[,1]==0,1]=NA
gamma_2_out_2[I_2_out_2[,2]==0,2]=NA
gamma_2_out_2[I_2_out_2[,3]==0,3]=NA
Tau_2_out_2=array(rtruncnorm(n=3*(N-n_0),a=d[D==1]-Tau_star_1,b=d[D==1],mean=d[D==1]-mu_tau_2_out_2,sd=sqrt(sigma_tau_2_2_out_2)),c(N-n_0,3))
Tau_2_out_2[I_2_out_2[,1]==0,1]=NA
Tau_2_out_2[I_2_out_2[,2]==0,2]=NA
Tau_2_out_2[I_2_out_2[,3]==0,3]=NA
sigma_2_2_out_2=3

#biomarker 3
mu_theta_3_out_2=rnorm(1,mean=mu_0_3,sd=sqrt(sigma_0_3_2))
sigma_theta_3_2_out_2=rinvgamma(1,shape=a_theta_3,scale=b_theta_3)
Theta_3_out_2=array(rnorm(N,mean=mu_theta_3_out_2,sd=sigma_theta_3_2_out_2),c(N,1))
I_3_out_2=array(rbinom(3*(N-n_0),size=1,prob=0.5),c(N-n_0,3))
mu_gamma_3_out_2=rnorm(1,mean=mu_1_3,sd=sqrt(sigma_1_3_2))
sigma_gamma_3_2_out_2=rinvgamma(1,shape=a_gamma_3,scale=b_gamma_3)
mu_tau_3_out_2=1
sigma_tau_3_2_out_2=(0.75)^2
gamma_3_out_2=array(exp(rnorm(3*(N-n_0),mean=mu_gamma_3_out_2,sd=sqrt(sigma_gamma_3_2_out_2))),c(N-n_0,3))
gamma_3_out_2[I_3_out_2[,1]==0,1]=NA
gamma_3_out_2[I_3_out_2[,2]==0,2]=NA
gamma_3_out_2[I_3_out_2[,3]==0,3]=NA
Tau_3_out_2=array(rtruncnorm(n=3*(N-n_0),a=d[D==1]-Tau_star_3,b=d[D==1],mean=d[D==1]-mu_tau_3_out_2,sd=sqrt(sigma_tau_3_2_out_2)),c(N-n_0,3))
Tau_3_out_2[I_3_out_2[,1]==0,1]=NA
Tau_3_out_2[I_3_out_2[,2]==0,2]=NA
Tau_3_out_2[I_3_out_2[,3]==0,3]=NA
sigma_3_2_out_2=1

#Update parameters
for(i in 1:S)
{
  #MU_THETA_1
  mu_theta_parameters=posterior_mu_theta(Theta=Theta_1_out_2[,i],sigma_theta_2=sigma_theta_1_2_out_2[i],mu_0=mu_0_1,sigma_0_2=sigma_0_1_2)
  mu_theta_1_out_2=c(mu_theta_1_out_2,rnorm(1,mean=mu_theta_parameters$mu_0_star,sd=sqrt(mu_theta_parameters$sigma_0_star_2)))
  
  #MU_THETA_2
  mu_theta_parameters=posterior_mu_theta(Theta=Theta_2_out_2[,i],sigma_theta_2=sigma_theta_2_2_out_2[i],mu_0=mu_0_2,sigma_0_2=sigma_0_2_2)
  mu_theta_2_out_2=c(mu_theta_2_out_2,rnorm(1,mean=mu_theta_parameters$mu_0_star,sd=sqrt(mu_theta_parameters$sigma_0_star_2)))
  
  #MU_THETA_3
  mu_theta_parameters=posterior_mu_theta(Theta=Theta_3_out_2[,i],sigma_theta_2=sigma_theta_3_2_out_2[i],mu_0=mu_0_3,sigma_0_2=sigma_0_3_2)
  mu_theta_3_out_2=c(mu_theta_3_out_2,rnorm(1,mean=mu_theta_parameters$mu_0_star,sd=sqrt(mu_theta_parameters$sigma_0_star_2)))
  
  #SIGMA_THETA_1
  sigma_theta_parameters=posterior_sigma_theta_2(Theta=Theta_1_out_2[,i],mu_theta=mu_theta_1_out_2[i+1],a_theta=a_theta_1,b_theta=b_theta_1)
  sigma_theta_1_2_out_2=c(sigma_theta_1_2_out_2,rinvgamma(1,shape=sigma_theta_parameters$a_theta_star,scale=sigma_theta_parameters$b_theta_star))
  
  #SIGMA_THETA_2
  sigma_theta_parameters=posterior_sigma_theta_2(Theta=Theta_2_out_2[,i],mu_theta=mu_theta_2_out_2[i+1],a_theta=a_theta_2,b_theta=b_theta_2)
  sigma_theta_2_2_out_2=c(sigma_theta_2_2_out_2,rinvgamma(1,shape=sigma_theta_parameters$a_theta_star,scale=sigma_theta_parameters$b_theta_star))
  
  #SIGMA_THETA_3
  sigma_theta_parameters=posterior_sigma_theta_2(Theta=Theta_3_out_2[,i],mu_theta=mu_theta_3_out_2[i+1],a_theta=a_theta_3,b_theta=b_theta_3)
  sigma_theta_3_2_out_2=c(sigma_theta_3_2_out_2,rinvgamma(1,shape=sigma_theta_parameters$a_theta_star,scale=sigma_theta_parameters$b_theta_star))
  
  #SIGMA_1
  sigma_parameters=posterior_sigma_2(Y=data$Y1,t=data$t,D=D,J=J,Theta=Theta_1_out_2[,i],I=I_1_out_2[,3*i],gamma=gamma_1_out_2[,3*i],Tau=Tau_1_out_2[,3*i])
  sigma_1_2_out_2=c(sigma_1_2_out_2,rinvgamma(1,shape=sigma_parameters$a_sigma,scale=sigma_parameters$b_sigma))
  
  #SIGMA_2
  sigma_parameters=posterior_sigma_2(Y=data$Y2,t=data$t,D=D,J=J,Theta=Theta_2_out_2[,i],I=I_2_out_2[,3*i],gamma=gamma_2_out_2[,3*i],Tau=Tau_2_out_2[,3*i])
  sigma_2_2_out_2=c(sigma_2_2_out_2,rinvgamma(1,shape=sigma_parameters$a_sigma,scale=sigma_parameters$b_sigma))
  
  #SIGMA_3
  sigma_parameters=posterior_sigma_2(Y=data$Y3,t=data$t,D=D,J=J,Theta=Theta_3_out_2[,i],I=I_3_out_2[,3*i],gamma=gamma_3_out_2[,3*i],Tau=Tau_3_out_2[,3*i])
  sigma_3_2_out_2=c(sigma_3_2_out_2,rinvgamma(1,shape=sigma_parameters$a_sigma,scale=sigma_parameters$b_sigma))
  
  #MU_I
  temp_I_matrix=cbind(I_1_out_2[,3*i],I_2_out_2[,3*i],I_3_out_2[,3*i])
  mu_I_output=update_mu_I(temp_I_matrix,mu_I_out_2[i],eta_I_out_2[i],delta_mu_I_set)
  mu_I_out_2=c(mu_I_out_2,mu_I_output$new_mu_I)
  
  #ETA_I
  temp_I_matrix=cbind(I_1_out_2[,3*i],I_2_out_2[,3*i],I_3_out_2[,3*i])
  eta_I_output=update_eta_I(temp_I_matrix,mu_I_out_2[i+1],eta_I_out_2[i],delta_eta_I_set)
  eta_I_out_2=c(eta_I_out_2,eta_I_output$new_eta_I)
  
  #MU_GAMMA_1
  mu_gamma_parameters=posterior_mu_gamma(I=I_1_out_2[,3*i],gamma=gamma_1_out_2[,3*i],sigma_gamma_2=sigma_gamma_1_2_out_2[i],mu_1=mu_1_1,sigma_1_2=sigma_1_1_2)
  mu_gamma_1_out_2=c(mu_gamma_1_out_2,rnorm(1,mean=mu_gamma_parameters$mu_1_star,sd=sqrt(mu_gamma_parameters$sigma_1_star_2)))
  
  #MU_GAMMA_2
  mu_gamma_parameters=posterior_mu_gamma(I=I_2_out_2[,3*i],gamma=gamma_2_out_2[,3*i],sigma_gamma_2=sigma_gamma_2_2_out_2[i],mu_1=mu_1_2,sigma_1_2=sigma_1_2_2)
  mu_gamma_2_out_2=c(mu_gamma_2_out_2,rnorm(1,mean=mu_gamma_parameters$mu_1_star,sd=sqrt(mu_gamma_parameters$sigma_1_star_2)))
  
  #MU_GAMMA_3
  mu_gamma_parameters=posterior_mu_gamma(I=I_3_out_2[,3*i],gamma=gamma_3_out_2[,3*i],sigma_gamma_2=sigma_gamma_3_2_out_2[i],mu_1=mu_1_3,sigma_1_2=sigma_1_3_2)
  mu_gamma_3_out_2=c(mu_gamma_3_out_2,rnorm(1,mean=mu_gamma_parameters$mu_1_star,sd=sqrt(mu_gamma_parameters$sigma_1_star_2)))
  
  #SIGMA_GAMMA_1
  sigma_gamma_parameters=posterior_sigma_gamma_2(I=I_1_out_2[,3*i],gamma=gamma_1_out_2[,3*i],mu_gamma=mu_gamma_1_out_2[i+1],a_gamma=a_gamma_1,b_gamma=b_gamma_1)
  sigma_gamma_1_2_out_2=c(sigma_gamma_1_2_out_2,rinvgamma(1,shape=sigma_gamma_parameters$a_gamma_star,scale=sigma_gamma_parameters$b_gamma_star))
  
  #SIGMA_GAMMA_2
  sigma_gamma_parameters=posterior_sigma_gamma_2(I=I_2_out_2[,3*i],gamma=gamma_2_out_2[,3*i],mu_gamma=mu_gamma_2_out_2[i+1],a_gamma=a_gamma_2,b_gamma=b_gamma_2)
  sigma_gamma_2_2_out_2=c(sigma_gamma_2_2_out_2,rinvgamma(1,shape=sigma_gamma_parameters$a_gamma_star,scale=sigma_gamma_parameters$b_gamma_star))
  
  #SIGMA_GAMMA_3
  sigma_gamma_parameters=posterior_sigma_gamma_2(I=I_3_out_2[,3*i],gamma=gamma_3_out_2[,3*i],mu_gamma=mu_gamma_3_out_2[i+1],a_gamma=a_gamma_3,b_gamma=b_gamma_3)
  sigma_gamma_3_2_out_2=c(sigma_gamma_3_2_out_2,rinvgamma(1,shape=sigma_gamma_parameters$a_gamma_star,scale=sigma_gamma_parameters$b_gamma_star))
  
  #MU_TAU_1
  mu_tau_output=update_mu_tau(Tau=Tau_1_out_2[,3*i],I=I_1_out_2[,3*i],d=d,D=D,mu_tau=mu_tau_1_out_2[i],sigma_tau_2=sigma_tau_1_2_out_2[i],delta_mu_tau=delta_mu_tau_set,mu_2=mu_2_1,sigma_2_2=sigma_2_1_2,Tau_star=Tau_star_1)
  mu_tau_1_out_2=c(mu_tau_1_out_2,mu_tau_output$new_mu_tau)
  
  #MU_TAU_2
  mu_tau_output=update_mu_tau(Tau=Tau_2_out_2[,3*i],I=I_2_out_2[,3*i],d=d,D=D,mu_tau=mu_tau_2_out_2[i],sigma_tau_2=sigma_tau_2_2_out_2[i],delta_mu_tau=delta_mu_tau_set,mu_2=mu_2_2,sigma_2_2=sigma_2_2_2,Tau_star=Tau_star_2)
  mu_tau_2_out_2=c(mu_tau_2_out_2,mu_tau_output$new_mu_tau)
  
  #MU_TAU_3
  mu_tau_output=update_mu_tau(Tau=Tau_3_out_2[,3*i],I=I_3_out_2[,3*i],d=d,D=D,mu_tau=mu_tau_3_out_2[i],sigma_tau_2=sigma_tau_3_2_out_2[i],delta_mu_tau=delta_mu_tau_set,mu_2=mu_2_3,sigma_2_2=sigma_2_3_2,Tau_star=Tau_star_3)
  mu_tau_3_out_2=c(mu_tau_3_out_2,mu_tau_output$new_mu_tau)
  
  #SIGMA_TAU_1
  sigma_tau_output=update_sigma_tau(Tau=Tau_1_out_2[,3*i],I=I_1_out_2[,3*i],d=d,D=D,mu_tau=mu_tau_1_out_2[i+1],sigma_tau_2=sigma_tau_1_2_out_2[i],delta_sigma_tau=delta_sigma_tau_set,a_tau=a_tau_1,b_tau=b_tau_1,Tau_star=Tau_star_1)
  sigma_tau_1_2_out_2=c(sigma_tau_1_2_out_2,sigma_tau_output$new_sigma_tau_2)
  
  #SIGMA_TAU_2
  sigma_tau_output=update_sigma_tau(Tau=Tau_2_out_2[,3*i],I=I_2_out_2[,3*i],d=d,D=D,mu_tau=mu_tau_2_out_2[i+1],sigma_tau_2=sigma_tau_2_2_out_2[i],delta_sigma_tau=delta_sigma_tau_set,a_tau=a_tau_2,b_tau=b_tau_2,Tau_star=Tau_star_2)
  sigma_tau_2_2_out_2=c(sigma_tau_2_2_out_2,sigma_tau_output$new_sigma_tau_2)
  
  #SIGMA_TAU_3
  sigma_tau_output=update_sigma_tau(Tau=Tau_3_out_2[,3*i],I=I_3_out_2[,3*i],d=d,D=D,mu_tau=mu_tau_3_out_2[i+1],sigma_tau_2=sigma_tau_3_2_out_2[i],delta_sigma_tau=delta_sigma_tau_set,a_tau=a_tau_3,b_tau=b_tau_3,Tau_star=Tau_star_3)
  sigma_tau_3_2_out_2=c(sigma_tau_3_2_out_2,sigma_tau_output$new_sigma_tau_2)
  
  Theta_1_new=rep(NA,N)
  Theta_2_new=rep(NA,N)
  Theta_3_new=rep(NA,N)
  
  I_1_new_1=rep(NA,N-n_0)
  I_1_new_2=rep(NA,N-n_0)
  I_1_new_3=rep(NA,N-n_0)
  
  I_2_new_1=rep(NA,N-n_0)
  I_2_new_2=rep(NA,N-n_0)
  I_2_new_3=rep(NA,N-n_0)
  
  I_3_new_1=rep(NA,N-n_0)
  I_3_new_2=rep(NA,N-n_0)
  I_3_new_3=rep(NA,N-n_0)
  
  gamma_1_new_1=rep(NA,N-n_0)
  gamma_1_new_2=rep(NA,N-n_0)
  gamma_1_new_3=rep(NA,N-n_0)
  
  gamma_2_new_1=rep(NA,N-n_0)
  gamma_2_new_2=rep(NA,N-n_0)
  gamma_2_new_3=rep(NA,N-n_0)
  
  gamma_3_new_1=rep(NA,N-n_0)
  gamma_3_new_2=rep(NA,N-n_0)
  gamma_3_new_3=rep(NA,N-n_0)
  
  Tau_1_new_1=rep(NA,N-n_0)
  Tau_1_new_2=rep(NA,N-n_0)
  Tau_1_new_3=rep(NA,N-n_0)
  
  Tau_2_new_1=rep(NA,N-n_0)
  Tau_2_new_2=rep(NA,N-n_0)
  Tau_2_new_3=rep(NA,N-n_0)
  
  Tau_3_new_1=rep(NA,N-n_0)
  Tau_3_new_2=rep(NA,N-n_0)
  Tau_3_new_3=rep(NA,N-n_0)
  
  for(j in 1:n_0)
  {
    #THETA_1_j
    theta_j_parameters=posterior_theta(Y_i=data$Y1[(data$ID==j)],t_i=data$t[(data$ID==j)],D_i=D[j],J_i=J[j],mu_theta=mu_theta_1_out_2[i+1],sigma_theta_2=sigma_theta_1_2_out_2[i+1],sigma_2=sigma_1_2_out_2[i+1],I_i=NA,gamma_i=NA,tau_i=NA)
    Theta_1_new[j]=rnorm(1,mean=theta_j_parameters$mu_theta_star,sd=sqrt(theta_j_parameters$sigma_theta_star_2))
    
    #THETA_2_j
    theta_j_parameters=posterior_theta(Y_i=data$Y2[(data$ID==j)],t_i=data$t[(data$ID==j)],D_i=D[j],J_i=J[j],mu_theta=mu_theta_2_out_2[i+1],sigma_theta_2=sigma_theta_2_2_out_2[i+1],sigma_2=sigma_2_2_out_2[i+1],I_i=NA,gamma_i=NA,tau_i=NA)
    Theta_2_new[j]=rnorm(1,mean=theta_j_parameters$mu_theta_star,sd=sqrt(theta_j_parameters$sigma_theta_star_2))
    
    #THETA_3_j
    theta_j_parameters=posterior_theta(Y_i=data$Y3[(data$ID==j)],t_i=data$t[(data$ID==j)],D_i=D[j],J_i=J[j],mu_theta=mu_theta_3_out_2[i+1],sigma_theta_2=sigma_theta_3_2_out_2[i+1],sigma_2=sigma_3_2_out_2[i+1],I_i=NA,gamma_i=NA,tau_i=NA)
    Theta_3_new[j]=rnorm(1,mean=theta_j_parameters$mu_theta_star,sd=sqrt(theta_j_parameters$sigma_theta_star_2))
  }
  for(j in (n_0+1):N)
  {
    #THETA_1_j
    theta_j_parameters=posterior_theta(Y_i=data$Y1[(data$ID==j)],t_i=data$t[(data$ID==j)],D_i=D[j],J_i=J[j],mu_theta=mu_theta_1_out_2[i+1],sigma_theta_2=sigma_theta_1_2_out_2[i+1],sigma_2=sigma_1_2_out_2[i+1],I_i=I_1_out_2[j-n_0,3*i],gamma_i=gamma_1_out_2[j-n_0,3*i],tau_i=Tau_1_out_2[j-n_0,3*i])
    Theta_1_new[j]=rnorm(1,mean=theta_j_parameters$mu_theta_star,sd=sqrt(theta_j_parameters$sigma_theta_star_2))
    
    #THETA_2_j
    theta_j_parameters=posterior_theta(Y_i=data$Y2[(data$ID==j)],t_i=data$t[(data$ID==j)],D_i=D[j],J_i=J[j],mu_theta=mu_theta_2_out_2[i+1],sigma_theta_2=sigma_theta_2_2_out_2[i+1],sigma_2=sigma_2_2_out_2[i+1],I_i=I_2_out_2[j-n_0,3*i],gamma_i=gamma_2_out_2[j-n_0,3*i],tau_i=Tau_2_out_2[j-n_0,3*i])
    Theta_2_new[j]=rnorm(1,mean=theta_j_parameters$mu_theta_star,sd=sqrt(theta_j_parameters$sigma_theta_star_2))
    
    #THETA_3_j
    theta_j_parameters=posterior_theta(Y_i=data$Y3[(data$ID==j)],t_i=data$t[(data$ID==j)],D_i=D[j],J_i=J[j],mu_theta=mu_theta_3_out_2[i+1],sigma_theta_2=sigma_theta_3_2_out_2[i+1],sigma_2=sigma_3_2_out_2[i+1],I_i=I_3_out_2[j-n_0,3*i],gamma_i=gamma_3_out_2[j-n_0,3*i],tau_i=Tau_3_out_2[j-n_0,3*i])
    Theta_3_new[j]=rnorm(1,mean=theta_j_parameters$mu_theta_star,sd=sqrt(theta_j_parameters$sigma_theta_star_2))
    
    #I_1
    F_temp=mu_I_out_2[i+1] + eta_I_out_2[i+1]*(I_2_out_2[j-n_0,3*i] + I_3_out_2[j-n_0,3*i])
    Pi_1_input=exp(F_temp)/(1+exp(F_temp))
    I_output=update_I(Y_i=data$Y1[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],d_i=d[j],mu_gamma=mu_gamma_1_out_2[i+1],sigma_gamma_2=sigma_gamma_1_2_out_2[i+1],mu_tau=mu_tau_1_out_2[i+1],sigma_tau_2=sigma_tau_1_2_out_2[i+1],Tau_star=Tau_star_1,sigma_2=sigma_1_2_out_2[i+1],Theta_i=Theta_1_new[j],I_i=I_1_out_2[j-n_0,3*i],gamma_i=gamma_1_out_2[j-n_0,3*i],tau_i=Tau_1_out_2[j-n_0,3*i],Pi=Pi_1_input)
    I_1_new_1[j-n_0]=I_output$new_I
    Tau_1_new_1[j-n_0]=I_output$new_tau
    gamma_1_new_1[j-n_0]=I_output$new_gamma
    
    #I_2
    F_temp=mu_I_out_2[i+1] + eta_I_out_2[i+1]*(I_1_new_1[j-n_0] + I_3_out_2[j-n_0,3*i])
    Pi_2_input=exp(F_temp)/(1+exp(F_temp))
    I_output=update_I(Y_i=data$Y2[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],d_i=d[j],mu_gamma=mu_gamma_2_out_2[i+1],sigma_gamma_2=sigma_gamma_2_2_out_2[i+1],mu_tau=mu_tau_2_out_2[i+1],sigma_tau_2=sigma_tau_2_2_out_2[i+1],Tau_star=Tau_star_2,sigma_2=sigma_2_2_out_2[i+1],Theta_i=Theta_2_new[j],I_i=I_2_out_2[j-n_0,3*i],gamma_i=gamma_2_out_2[j-n_0,3*i],tau_i=Tau_2_out_2[j-n_0,3*i],Pi=Pi_2_input)
    I_2_new_1[j-n_0]=I_output$new_I
    Tau_2_new_1[j-n_0]=I_output$new_tau
    gamma_2_new_1[j-n_0]=I_output$new_gamma
    
    #I_3
    F_temp=mu_I_out_2[i+1] + eta_I_out_2[i+1]*(I_1_new_1[j-n_0] + I_2_new_1[j-n_0])
    Pi_3_input=exp(F_temp)/(1+exp(F_temp))
    I_output=update_I(Y_i=data$Y3[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],d_i=d[j],mu_gamma=mu_gamma_3_out_2[i+1],sigma_gamma_2=sigma_gamma_3_2_out_2[i+1],mu_tau=mu_tau_3_out_2[i+1],sigma_tau_2=sigma_tau_3_2_out_2[i+1],Tau_star=Tau_star_3,sigma_2=sigma_3_2_out_2[i+1],Theta_i=Theta_3_new[j],I_i=I_3_out_2[j-n_0,3*i],gamma_i=gamma_3_out_2[j-n_0,3*i],tau_i=Tau_3_out_2[j-n_0,3*i],Pi=Pi_3_input)
    I_3_new_1[j-n_0]=I_output$new_I
    Tau_3_new_1[j-n_0]=I_output$new_tau
    gamma_3_new_1[j-n_0]=I_output$new_gamma
    
    #GAMMA_1
    gamma_output=update_gamma(Y_i=data$Y1[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],mu_gamma=mu_gamma_1_out_2[i+1],sigma_gamma_2=sigma_gamma_1_2_out_2[i+1],sigma_2=sigma_1_2_out_2[i+1],Theta_i=Theta_1_new[j],I_i=I_1_new_1[j-n_0],gamma_i=gamma_1_new_1[j-n_0],tau_i=Tau_1_new_1[j-n_0],delta_gamma=delta_gamma_set)
    gamma_1_new_2[j-n_0]=gamma_output$new_gamma
    I_1_new_2[j-n_0]=gamma_output$new_I
    Tau_1_new_2[j-n_0]=gamma_output$new_tau
    
    #GAMMA_2
    gamma_output=update_gamma(Y_i=data$Y2[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],mu_gamma=mu_gamma_2_out_2[i+1],sigma_gamma_2=sigma_gamma_2_2_out_2[i+1],sigma_2=sigma_2_2_out_2[i+1],Theta_i=Theta_2_new[j],I_i=I_2_new_1[j-n_0],gamma_i=gamma_2_new_1[j-n_0],tau_i=Tau_2_new_1[j-n_0],delta_gamma=delta_gamma_set)
    gamma_2_new_2[j-n_0]=gamma_output$new_gamma
    I_2_new_2[j-n_0]=gamma_output$new_I
    Tau_2_new_2[j-n_0]=gamma_output$new_tau
    
    #GAMMA_3
    gamma_output=update_gamma(Y_i=data$Y3[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],mu_gamma=mu_gamma_3_out_2[i+1],sigma_gamma_2=sigma_gamma_3_2_out_2[i+1],sigma_2=sigma_3_2_out_2[i+1],Theta_i=Theta_3_new[j],I_i=I_3_new_1[j-n_0],gamma_i=gamma_3_new_1[j-n_0],tau_i=Tau_3_new_1[j-n_0],delta_gamma=delta_gamma_set)
    gamma_3_new_2[j-n_0]=gamma_output$new_gamma
    I_3_new_2[j-n_0]=gamma_output$new_I
    Tau_3_new_2[j-n_0]=gamma_output$new_tau
    
    #TAU_1
    tau_output=update_tau(Y_i=data$Y1[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],d_i=d[j],sigma_2=sigma_1_2_out_2[i+1],mu_tau=mu_tau_1_out_2[i+1],sigma_tau_2=sigma_tau_1_2_out_2[i+1],Theta_i=Theta_1_new[j],I_i=I_1_new_2[j-n_0],gamma_i=gamma_1_new_2[j-n_0],tau_i=Tau_1_new_2[j-n_0],delta_tau=delta_tau_set,Tau_star=Tau_star_1)
    I_1_new_3[j-n_0]=tau_output$new_I
    Tau_1_new_3[j-n_0]=tau_output$new_tau
    gamma_1_new_3[j-n_0]=tau_output$new_gamma
    
    #TAU_2
    tau_output=update_tau(Y_i=data$Y2[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],d_i=d[j],sigma_2=sigma_2_2_out_2[i+1],mu_tau=mu_tau_2_out_2[i+1],sigma_tau_2=sigma_tau_2_2_out_2[i+1],Theta_i=Theta_2_new[j],I_i=I_2_new_2[j-n_0],gamma_i=gamma_2_new_2[j-n_0],tau_i=Tau_2_new_2[j-n_0],delta_tau=delta_tau_set,Tau_star=Tau_star_2)
    I_2_new_3[j-n_0]=tau_output$new_I
    Tau_2_new_3[j-n_0]=tau_output$new_tau
    gamma_2_new_3[j-n_0]=tau_output$new_gamma
    
    #TAU_1
    tau_output=update_tau(Y_i=data$Y3[(data$ID==j)],t_i=data$t[(data$ID==j)],J_i=J[j],d_i=d[j],sigma_2=sigma_3_2_out_2[i+1],mu_tau=mu_tau_3_out_2[i+1],sigma_tau_2=sigma_tau_3_2_out_2[i+1],Theta_i=Theta_3_new[j],I_i=I_3_new_2[j-n_0],gamma_i=gamma_3_new_2[j-n_0],tau_i=Tau_3_new_2[j-n_0],delta_tau=delta_tau_set,Tau_star=Tau_star_3)
    I_3_new_3[j-n_0]=tau_output$new_I
    Tau_3_new_3[j-n_0]=tau_output$new_tau
    gamma_3_new_3[j-n_0]=tau_output$new_gamma
    
  }
  Theta_1_out_2=cbind(Theta_1_out_2,Theta_1_new)
  Theta_2_out_2=cbind(Theta_2_out_2,Theta_2_new)
  Theta_3_out_2=cbind(Theta_3_out_2,Theta_3_new)
  
  I_1_out_2=cbind(I_1_out_2,I_1_new_1,I_1_new_2,I_1_new_3)
  I_2_out_2=cbind(I_2_out_2,I_2_new_1,I_2_new_2,I_2_new_3)
  I_3_out_2=cbind(I_3_out_2,I_3_new_1,I_3_new_2,I_3_new_3)
  
  gamma_1_out_2=cbind(gamma_1_out_2,gamma_1_new_1,gamma_1_new_2,gamma_1_new_3)
  gamma_2_out_2=cbind(gamma_2_out_2,gamma_2_new_1,gamma_2_new_2,gamma_2_new_3)
  gamma_3_out_2=cbind(gamma_3_out_2,gamma_3_new_1,gamma_3_new_2,gamma_3_new_3)
  
  Tau_1_out_2=cbind(Tau_1_out_2,Tau_1_new_1,Tau_1_new_2,Tau_1_new_3)
  Tau_2_out_2=cbind(Tau_2_out_2,Tau_2_new_1,Tau_2_new_2,Tau_2_new_3)
  Tau_3_out_2=cbind(Tau_3_out_2,Tau_3_new_1,Tau_3_new_2,Tau_3_new_3)
  
  print(i)
}

########################################################################################################################################################################################################################
#Post-processing results from training data
########################################################################################################################################################################################################################
#TRACE PLOTS
trace_plot=function(y1,y2,title)
{
  ymin=min(y1,y2)
  ymax=max(y1,y2)
  x=seq(from=1,to=length(y1),by=1)
  plot(x,y1,type="l",xlab="Iteration Number",ylab=title,col="red",ylim=c(ymin,ymax))
  x=seq(from=1,to=length(y2),by=1)
  lines(x,y2,lty=1,col="blue")
}

trace_plot(mu_theta_1_out,mu_theta_1_out_2,"MU_THETA_1")
trace_plot(mu_theta_2_out,mu_theta_2_out_2,"MU_THETA_2")
trace_plot(mu_theta_3_out,mu_theta_3_out_2,"MU_THETA_3")

trace_plot(sigma_theta_1_2_out,sigma_theta_1_2_out_2,"SIGMA_THETA_1_2")
trace_plot(sigma_theta_2_2_out,sigma_theta_2_2_out_2,"SIGMA_THETA_2_2")
trace_plot(sigma_theta_3_2_out,sigma_theta_3_2_out_2,"SIGMA_THETA_3_2")

trace_plot(sigma_1_2_out,sigma_1_2_out_2,"SIGMA_1_2")
trace_plot(sigma_2_2_out,sigma_2_2_out_2,"SIGMA_2_2")
trace_plot(sigma_3_2_out,sigma_3_2_out_2,"SIGMA_3_2")

trace_plot(mu_I_out,mu_I_out_2,"MU_I")
trace_plot(eta_I_out,eta_I_out_2,"ETA_I")

trace_plot(mu_gamma_1_out,mu_gamma_1_out_2,"MU_GAMMA_1")
trace_plot(mu_gamma_2_out,mu_gamma_2_out_2,"MU_GAMMA_2")
trace_plot(mu_gamma_3_out,mu_gamma_3_out_2,"MU_GAMMA_3")

trace_plot(sigma_gamma_1_2_out,sigma_gamma_1_2_out_2,"SIGMA_GAMMA_1_2")
trace_plot(sigma_gamma_2_2_out,sigma_gamma_2_2_out_2,"SIGMA_GAMMA_2_2")
trace_plot(sigma_gamma_3_2_out,sigma_gamma_3_2_out_2,"SIGMA_GAMMA_3_2")

trace_plot(mu_tau_1_out,mu_tau_1_out_2,"MU_TAU_1")
trace_plot(mu_tau_2_out,mu_tau_2_out_2,"MU_TAU_2")
trace_plot(mu_tau_3_out,mu_tau_3_out_2,"MU_TAU_3")

trace_plot(sigma_tau_1_2_out,sigma_tau_1_2_out_2,"SIGMA_TAU_1_2")
trace_plot(sigma_tau_2_2_out,sigma_tau_2_2_out_2,"SIGMA_TAU_2_2")
trace_plot(sigma_tau_3_2_out,sigma_tau_3_2_out_2,"SIGMA_TAU_3_2")

n_thin=10 #Thin to reduce autocorrelation
B=2000 #burn-in interations
B_star=3*B

mu_theta_1_out_thin=mu_theta_1_out[seq(from=B+1,to=S+1,by=n_thin)]
mu_theta_1_out_2_thin=mu_theta_1_out_2[seq(from=B+1,to=S+1,by=n_thin)]
mu_theta_2_out_thin=mu_theta_2_out[seq(from=B+1,to=S+1,by=n_thin)]
mu_theta_2_out_2_thin=mu_theta_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]
mu_theta_3_out_thin=mu_theta_3_out[seq(from=B+1,to=S+1,by=n_thin)]
mu_theta_3_out_2_thin=mu_theta_3_out_2[seq(from=B+1,to=S+1,by=n_thin)]

sigma_theta_1_2_out_thin=sigma_theta_1_2_out[seq(from=B+1,to=S+1,by=n_thin)]
sigma_theta_1_2_out_2_thin=sigma_theta_1_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]
sigma_theta_2_2_out_thin=sigma_theta_2_2_out[seq(from=B+1,to=S+1,by=n_thin)]
sigma_theta_2_2_out_2_thin=sigma_theta_2_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]
sigma_theta_3_2_out_thin=sigma_theta_3_2_out[seq(from=B+1,to=S+1,by=n_thin)]
sigma_theta_3_2_out_2_thin=sigma_theta_3_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]

sigma_1_2_out_thin=sigma_1_2_out[seq(from=B+1,to=S+1,by=n_thin)]
sigma_1_2_out_2_thin=sigma_1_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]
sigma_2_2_out_thin=sigma_2_2_out[seq(from=B+1,to=S+1,by=n_thin)]
sigma_2_2_out_2_thin=sigma_2_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]
sigma_3_2_out_thin=sigma_3_2_out[seq(from=B+1,to=S+1,by=n_thin)]
sigma_3_2_out_2_thin=sigma_3_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]

mu_I_out_thin=mu_I_out[seq(from=B+1,to=S+1,by=n_thin)]
mu_I_out_2_thin=mu_I_out_2[seq(from=B+1,to=S+1,by=n_thin)]

eta_I_out_thin=eta_I_out[seq(from=B+1,to=S+1,by=n_thin)]
eta_I_out_2_thin=eta_I_out_2[seq(from=B+1,to=S+1,by=n_thin)]

mu_gamma_1_out_thin=mu_gamma_1_out[seq(from=B+1,to=S+1,by=n_thin)]
mu_gamma_1_out_2_thin=mu_gamma_1_out_2[seq(from=B+1,to=S+1,by=n_thin)]
mu_gamma_2_out_thin=mu_gamma_2_out[seq(from=B+1,to=S+1,by=n_thin)]
mu_gamma_2_out_2_thin=mu_gamma_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]
mu_gamma_3_out_thin=mu_gamma_3_out[seq(from=B+1,to=S+1,by=n_thin)]
mu_gamma_3_out_2_thin=mu_gamma_3_out_2[seq(from=B+1,to=S+1,by=n_thin)]

sigma_gamma_1_2_out_thin=sigma_gamma_1_2_out[seq(from=B+1,to=S+1,by=n_thin)]
sigma_gamma_1_2_out_2_thin=sigma_gamma_1_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]
sigma_gamma_2_2_out_thin=sigma_gamma_2_2_out[seq(from=B+1,to=S+1,by=n_thin)]
sigma_gamma_2_2_out_2_thin=sigma_gamma_2_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]
sigma_gamma_3_2_out_thin=sigma_gamma_3_2_out[seq(from=B+1,to=S+1,by=n_thin)]
sigma_gamma_3_2_out_2_thin=sigma_gamma_3_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]

mu_tau_1_out_thin=mu_tau_1_out[seq(from=B+1,to=S+1,by=n_thin)]
mu_tau_1_out_2_thin=mu_tau_1_out_2[seq(from=B+1,to=S+1,by=n_thin)]
mu_tau_2_out_thin=mu_tau_2_out[seq(from=B+1,to=S+1,by=n_thin)]
mu_tau_2_out_2_thin=mu_tau_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]
mu_tau_3_out_thin=mu_tau_3_out[seq(from=B+1,to=S+1,by=n_thin)]
mu_tau_3_out_2_thin=mu_tau_3_out_2[seq(from=B+1,to=S+1,by=n_thin)]

sigma_tau_1_2_out_thin=sigma_tau_1_2_out[seq(from=B+1,to=S+1,by=n_thin)]
sigma_tau_1_2_out_2_thin=sigma_tau_1_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]
sigma_tau_2_2_out_thin=sigma_tau_2_2_out[seq(from=B+1,to=S+1,by=n_thin)]
sigma_tau_2_2_out_2_thin=sigma_tau_2_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]
sigma_tau_3_2_out_thin=sigma_tau_3_2_out[seq(from=B+1,to=S+1,by=n_thin)]
sigma_tau_3_2_out_2_thin=sigma_tau_3_2_out_2[seq(from=B+1,to=S+1,by=n_thin)]

#posterior_chain=an array with a column for each chain with burn-in removed
gelman_rubin_statistic=function(posterior_chain)
{
  m_star=ncol(posterior_chain)
  n_star=nrow(posterior_chain)
  m=m_star*2
  n=floor(n_star/2)
  posterior_chain_2=cbind(posterior_chain[1:n,],posterior_chain[(n+1):(2*n),])
  posterior_chain_2_mean=apply(posterior_chain_2,2,mean)
  B_GR=n*var(posterior_chain_2_mean)
  posterior_chain_2_var=apply(posterior_chain_2,2,var)
  W_GR=mean(posterior_chain_2_var)
  R_hat=sqrt(((n-1)*W_GR/n + B_GR/n)/W_GR)
  R_hat
}

summary_stats=NULL
#MU_THETA_1
title="MU_THETA_1"
x1=mu_theta_1_out_thin
x2=mu_theta_1_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#MU_THETA_2
title="MU_THETA_2"
x1=mu_theta_2_out_thin
x2=mu_theta_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#MU_THETA_3
title="MU_THETA_3"
x1=mu_theta_3_out_thin
x2=mu_theta_3_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#SIGMA_THETA_1
title="SIGMA_THETA_1_2"
x1=sigma_theta_1_2_out_thin
x2=sigma_theta_1_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#SIGMA_THETA_2
title="SIGMA_THETA_2_2"
x1=sigma_theta_2_2_out_thin
x2=sigma_theta_2_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#SIGMA_THETA_3
title="SIGMA_THETA_3_2"
x1=sigma_theta_3_2_out_thin
x2=sigma_theta_3_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#SIGMA_1
title="SIGMA_1_2"
x1=sigma_1_2_out_thin
x2=sigma_1_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#SIGMA_2
title="SIGMA_2_2"
x1=sigma_2_2_out_thin
x2=sigma_2_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#SIGMA_3
title="SIGMA_3_2"
x1=sigma_3_2_out_thin
x2=sigma_3_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#MU_I
title="MU_I"
x1=mu_I_out_thin
x2=mu_I_out_2_thin
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#ETA_I
title="ETA_I"
x1=eta_I_out_thin
x2=eta_I_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#MU_GAMMA_1
title="MU_GAMMA_1"
x1=mu_gamma_1_out_thin
x2=mu_gamma_1_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#MU_GAMMA_2
title="MU_GAMMA_2"
x1=mu_gamma_2_out_thin
x2=mu_gamma_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#MU_GAMMA_3
title="MU_GAMMA_3"
x1=mu_gamma_3_out_thin
x2=mu_gamma_3_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#SIGMA_GAMMA_1
title="SIGMA_GAMMA_1_2"
x1=sigma_gamma_1_2_out_thin
x2=sigma_gamma_1_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#SIGMA_GAMMA_2
title="SIGMA_GAMMA_2_2"
x1=sigma_gamma_2_2_out_thin
x2=sigma_gamma_2_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#SIGMA_GAMMA_3
title="SIGMA_GAMMA_3_2"
x1=sigma_gamma_3_2_out_thin
x2=sigma_gamma_3_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#MU_TAU_1
title="MU_TAU_1"
x1=mu_tau_1_out_thin
x2=mu_tau_1_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#MU_TAU_2
title="MU_TAU_2"
x1=mu_tau_2_out_thin
x2=mu_tau_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#MU_TAU_3
title="MU_TAU_3"
x1=mu_tau_3_out_thin
x2=mu_tau_3_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#SIGMA_TAU_1
title="SIGMA_TAU_1_2"
x1=sigma_tau_1_2_out_thin
x2=sigma_tau_1_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#SIGMA_TAU_2
title="SIGMA_TAU_2_2"
x1=sigma_tau_2_2_out_thin
x2=sigma_tau_2_2_out_2_thin
acf(x1,lag.max=100,main=paste("Chain 1",title,sep="-"))
acf(x2,lag.max=100,main=paste("Chain 2",title,sep="-"))
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

#SIGMA_TAU_3
title="SIGMA_TAU_3_2"
x1=sigma_tau_3_2_out_thin
x2=sigma_tau_3_2_out_2_thin
x=c(x1,x2)
summary_stats=rbind(summary_stats,cbind(Parameter=title,Mean=mean(x),Median=median(x),Quantile_025=as.numeric(quantile(x,probs=0.025)),Quantile_975=as.numeric(quantile(x,probs=0.975)),Min=min(x),Max=max(x),R_hat=gelman_rubin_statistic(cbind(x1,x2))))

## Output of summary_stats file:
#       Parameter         Mean                 Median               Quantile_025         Quantile_975         Min                   Max                 R_hat              
# [1,] "MU_THETA_1"      "2.42447964307287"   "2.42376713291975"   "2.33357590043021"   "2.51572210935502"   "2.2661551119047"     "2.58122274790118"  "0.999027055490734"
# [2,] "MU_THETA_2"      "3.08818602111016"   "3.08756319224239"   "2.98444315383398"   "3.19562932139037"   "2.91691392076204"    "3.27407974692887"  "0.999595209305632"
# [3,] "MU_THETA_3"      "2.69062812555171"   "2.69102549704294"   "2.59304611599708"   "2.7860173288433"    "2.53339536065743"    "2.86983558456411"  "1.00085689942108" 
# [4,] "SIGMA_THETA_1_2" "0.819127686520992"  "0.815872182627638"  "0.707110673642422"  "0.94856595515212"   "0.653450775133136"   "1.08867587606734"  "0.999766175425271"
# [5,] "SIGMA_THETA_2_2" "0.801820620837849"  "0.798369578287224"  "0.661535601276949"  "0.965215546728426"  "0.558142149538541"   "1.08466042376793"  "1.00001390421679" 
# [6,] "SIGMA_THETA_3_2" "0.741516111105822"  "0.735341493994697"  "0.619304239356033"  "0.885206051026522"  "0.569354303052035"   "0.979353971236852" "1.00174215438518" 
# [7,] "SIGMA_1_2"       "0.230468569677368"  "0.230474209285454"  "0.215738799555705"  "0.246133584092399"  "0.203896647181668"   "0.258904292234748" "1.00134850930238" 
# [8,] "SIGMA_2_2"       "1.41106083897877"   "1.41071691931127"   "1.32037337428281"   "1.50230536876237"   "1.27042498168633"    "1.5878378501068"   "0.999029924206092"
# [9,] "SIGMA_3_2"       "0.820333991058462"  "0.820130967790534"  "0.769207008168376"  "0.875227384750189"  "0.733435500579164"   "0.910179328192249" "0.999147922342526"
# [10,] "MU_I"            "-0.16723249184188"  "-0.169955568768551" "-0.516920273175212" "0.203133455394032"  "-0.763199134550707"  "0.424871988456401" "0.999796623221603"
# [11,] "ETA_I"           "0.0957983880044332" "0.091765083805436"  "0.0321361865798529" "0.186703919161576"  "0.0164688073973981"  "0.244996683375926" "1.00010350073318" 
# [12,] "MU_GAMMA_1"      "2.39664637317449"   "2.39800748611611"   "2.0364734647599"    "2.76169456583781"   "1.71233128615782"    "3.0757829308187"   "1.02832948462635" 
# [13,] "MU_GAMMA_2"      "1.99312474865304"   "1.99596561161928"   "1.78555362812693"   "2.1975799613881"    "1.65701871445047"    "2.44519176282092"  "1.02273491538543" 
# [14,] "MU_GAMMA_3"      "1.82630767027032"   "1.81950524785114"   "1.50077908245071"   "2.23050557127802"   "1.37391557032121"    "2.48162830593404"  "1.05660166746905" 
# [15,] "SIGMA_GAMMA_1_2" "0.654679331511895"  "0.602116393248257"  "0.299059832685417"  "1.25875098860018"   "0.220265138126012"   "2.79871504378628"  "1.01119791561073" 
# [16,] "SIGMA_GAMMA_2_2" "0.0350995015308889" "0.0286355097737781" "0.0124566181937708" "0.0930956019012701" "0.00699129846518533" "0.386879635507194" "1.01305766651943" 
# [17,] "SIGMA_GAMMA_3_2" "0.0433543243763372" "0.0320364581352531" "0.0117481384421143" "0.148017465711593"  "0.00596956901239703" "0.503487893049127" "1.00939443382486" 
# [18,] "MU_TAU_1"        "1.00141493722285"   "0.998878957764709"  "0.524289504443929"  "1.48523014189714"   "0.146787272783898"   "1.92596958681683"  "1.02940532791427" 
# [19,] "MU_TAU_2"        "0.822879813996256"  "0.825201196527897"  "0.306610438377672"  "1.32347850156948"   "-0.10959252988393"   "1.69211833666668"  "1.00408250026214" 
# [20,] "MU_TAU_3"        "0.6169322159136"    "0.61315421758716"   "0.154239561907985"  "1.12684475929627"   "-0.139714728041973"  "1.6917253245588"   "1.0108945350872"  
# [21,] "SIGMA_TAU_1_2"   "0.710514246320026"  "0.665343136902194"  "0.382521243772831"  "1.33017059513714"   "0.255406161443289"   "2.58526165785156"  "1.00166667719488" 
# [22,] "SIGMA_TAU_2_2"   "0.687447140284648"  "0.64357983119657"   "0.363117437214555"  "1.28502829751126"   "0.23294276204632"    "2.30545128712923"  "0.999734517090155"
# [23,] "SIGMA_TAU_3_2"   "0.555209676563166"  "0.522430001626923"  "0.308231296931259"  "1.05556836898524"   "0.205458842232671"   "1.69237532161211"  "1.00206315321289" 


#Summary of Gelman-Rubin statistics
sum(summary_stats[,8]>1.01)
summary_stats[summary_stats[,8]>1.01,c(1,8)]
#All less than 1.06

#Extract samples from posterior distributions of hyper parameters to use in implementation of screening in the validation cohort
sigma_1_posterior=c(sigma_1_2_out_thin,sigma_1_2_out_2_thin)
sigma_2_posterior=c(sigma_2_2_out_thin,sigma_2_2_out_2_thin)
sigma_3_posterior=c(sigma_3_2_out_thin,sigma_3_2_out_2_thin)
mu_theta_1_posterior=c(mu_theta_1_out_thin,mu_theta_1_out_2_thin)
mu_theta_2_posterior=c(mu_theta_2_out_thin,mu_theta_2_out_2_thin)
mu_theta_3_posterior=c(mu_theta_3_out_thin,mu_theta_3_out_2_thin)
sigma_theta_1_2_posterior=c(sigma_theta_1_2_out_thin,sigma_theta_1_2_out_2_thin)
sigma_theta_2_2_posterior=c(sigma_theta_2_2_out_thin,sigma_theta_2_2_out_2_thin)
sigma_theta_3_2_posterior=c(sigma_theta_3_2_out_thin,sigma_theta_3_2_out_2_thin)
mu_I_posterior=c(mu_I_out_thin,mu_I_out_2_thin)
eta_I_posterior=c(eta_I_out_thin,eta_I_out_2_thin)
mu_gamma_1_posterior=c(mu_gamma_1_out_thin,mu_gamma_1_out_2_thin)
mu_gamma_2_posterior=c(mu_gamma_2_out_thin,mu_gamma_2_out_2_thin)
mu_gamma_3_posterior=c(mu_gamma_3_out_thin,mu_gamma_3_out_2_thin)
sigma_gamma_1_2_posterior=c(sigma_gamma_1_2_out_thin,sigma_gamma_1_2_out_2_thin)
sigma_gamma_2_2_posterior=c(sigma_gamma_2_2_out_thin,sigma_gamma_2_2_out_2_thin)
sigma_gamma_3_2_posterior=c(sigma_gamma_3_2_out_thin,sigma_gamma_3_2_out_2_thin)
mu_tau_1_posterior=c(mu_tau_1_out_thin,mu_tau_1_out_2_thin)
mu_tau_2_posterior=c(mu_tau_2_out_thin,mu_tau_2_out_2_thin)
mu_tau_3_posterior=c(mu_tau_3_out_thin,mu_tau_3_out_2_thin)
sigma_tau_1_2_posterior=c(sigma_tau_1_2_out_thin,sigma_tau_1_2_out_2_thin)
sigma_tau_2_2_posterior=c(sigma_tau_2_2_out_thin,sigma_tau_2_2_out_2_thin)
sigma_tau_3_2_posterior=c(sigma_tau_3_2_out_thin,sigma_tau_3_2_out_2_thin)
S=length(mu_I_posterior)

########################################################################################################################################################################################################################
#Implement screening in the validation data (pg 11 Supplementary Materials)
########################################################################################################################################################################################################################
#Get sample size 
N=length(unique(data_validation$ID)) #total number of patients
n_0=length(unique(data_validation$ID[(data_validation$D==0)])) #number of control patients

#Extract vectors
d=data_validation$d[(data_validation$obs_number==1)]
D=data_validation$D[(data_validation$obs_number==1)]
J=rle(data_validation$ID)$lengths
subject_ID=unique(data_validation$ID)

#observed HCC diagnosis times in training cohort
d_HCC=sort(unique(data$d[data$D==1]))

p_Y1_noHCC=rep(NA,length(data_validation$ID))
p_Y2_noHCC=rep(NA,length(data_validation$ID))
p_Y3_noHCC=rep(NA,length(data_validation$ID))
p_Y1_HCC=rep(NA,length(data_validation$ID))
p_Y2_HCC=rep(NA,length(data_validation$ID))
p_Y3_HCC=rep(NA,length(data_validation$ID))

l=1
for(i in 1:N)
{
  subject_data=subset(data_validation,data_validation$ID==subject_ID[i]) 
  
  for(j in 1:J[i])
  {    
    temp_Y1=rep(subject_data$Y1[1:j],rep(S,j))
    temp_Y2=rep(subject_data$Y2[1:j],rep(S,j))
    temp_Y3=rep(subject_data$Y3[1:j],rep(S,j))
    temp_t=rep(subject_data$t[1:j],rep(S,j))
    
    #Calculate Pr(Y1|No HCC) for each patient at each screening time
    theta_1_posterior=rnorm(n=S,mean=mu_theta_1_posterior,sd=sqrt(sigma_theta_1_2_posterior))
    
    prob_Y1_noHCC=dnorm(temp_Y1,mean=theta_1_posterior,sd=sqrt(sigma_1_posterior))
    prob_Y1_noHCC_1=array(prob_Y1_noHCC,c(S,j))
    prob_Y1_noHCC_2=apply(prob_Y1_noHCC_1,1,prod,na.rm=TRUE)
    p_Y1_noHCC[l]=mean(prob_Y1_noHCC_2)
    
    #Calculate Pr(Y2|No HCC) for each patient at each screening time
    theta_2_posterior=rnorm(n=S,mean=mu_theta_2_posterior,sd=sqrt(sigma_theta_2_2_posterior))
    
    prob_Y2_noHCC=dnorm(temp_Y2,mean=theta_2_posterior,sd=sqrt(sigma_2_posterior))
    prob_Y2_noHCC_1=array(prob_Y2_noHCC,c(S,j))
    prob_Y2_noHCC_2=apply(prob_Y2_noHCC_1,1,prod,na.rm=TRUE)
    p_Y2_noHCC[l]=mean(prob_Y2_noHCC_2)
    
    #Calculate Pr(Y3|No HCC) for each patient at each screening time
    theta_3_posterior=rnorm(n=S,mean=mu_theta_3_posterior,sd=sqrt(sigma_theta_3_2_posterior))
    
    prob_Y3_noHCC=dnorm(temp_Y3,mean=theta_3_posterior,sd=sqrt(sigma_3_posterior))
    prob_Y3_noHCC_1=array(prob_Y3_noHCC,c(S,j))
    prob_Y3_noHCC_2=apply(prob_Y3_noHCC_1,1,prod,na.rm=TRUE)
    p_Y3_noHCC[l]=mean(prob_Y3_noHCC_2)
    
    #Get draws of vector c(I_1,I_2)
    temp_p_1=get_p1(mu_I=mu_I_posterior,eta_I=eta_I_posterior)
    I_1_posterior=rbinom(n=S,size=1,prob=temp_p_1)
    temp_p_2=get_p2_I_1(I_1=I_1_posterior,mu_I=mu_I_posterior,eta_I=eta_I_posterior)
    I_2_posterior=rbinom(n=S,size=1,prob=temp_p_2)
    temp_p_3=get_p3_I_1_2(I_1=I_1_posterior,I_2=I_2_posterior,mu_I=mu_I_posterior,eta_I=eta_I_posterior)
    I_3_posterior=rbinom(n=S,size=1,prob=temp_p_3)
    
    ############ Draw d from empirical distribution ############
    u_temp=runif(S)
    d_draws=quantile(d_HCC,u_temp,type=1)
    
    #Calculate Pr(Y1|HCC) for each patient at each screening time
    theta_1_posterior=rnorm(n=S,mean=mu_theta_1_posterior,sd=sqrt(sigma_theta_1_2_posterior))
    gamma_1_posterior=exp(rnorm(n=S,mean=mu_gamma_1_posterior,sd=sqrt(sigma_gamma_1_2_posterior)))
    Tau_1_posterior=rtruncnorm(n=S,a=(d_draws-Tau_star_1),b=d_draws,mean=d_draws-mu_tau_1_posterior,sd=sqrt(sigma_tau_1_2_posterior))    
    mean_Y1=rep(theta_1_posterior,j) + rep(I_1_posterior,j)*rep(gamma_1_posterior,j)*(temp_t-Tau_1_posterior)*as.numeric(temp_t>Tau_1_posterior)
    
    prob_Y1_HCC=dnorm(temp_Y1,mean=mean_Y1,sd=sqrt(sigma_1_posterior))
    prob_Y1_HCC_1=array(prob_Y1_HCC,c(S,j))
    prob_Y1_HCC_2=apply(prob_Y1_HCC_1,1,prod,na.rm=TRUE)
    p_Y1_HCC[l]=mean(prob_Y1_HCC_2)
    
    #Calculate Pr(Y2|HCC) for each patient at each screening time
    theta_2_posterior=rnorm(n=S,mean=mu_theta_2_posterior,sd=sqrt(sigma_theta_2_2_posterior))
    gamma_2_posterior=exp(rnorm(n=S,mean=mu_gamma_2_posterior,sd=sqrt(sigma_gamma_2_2_posterior)))
    Tau_2_posterior=rtruncnorm(n=S,a=(d_draws-Tau_star_2),b=d_draws,mean=d_draws-mu_tau_2_posterior,sd=sqrt(sigma_tau_2_2_posterior))
    mean_Y2=rep(theta_2_posterior,j) + rep(I_2_posterior,j)*rep(gamma_2_posterior,j)*(temp_t-Tau_2_posterior)*as.numeric(temp_t>Tau_2_posterior)
    
    prob_Y2_HCC=dnorm(temp_Y2,mean=mean_Y2,sd=sqrt(sigma_2_posterior))
    prob_Y2_HCC_1=array(prob_Y2_HCC,c(S,j))
    prob_Y2_HCC_2=apply(prob_Y2_HCC_1,1,prod,na.rm=TRUE)
    p_Y2_HCC[l]=mean(prob_Y2_HCC_2)
    
    #Calculate Pr(Y3|HCC) for each patient at each screening time
    theta_3_posterior=rnorm(n=S,mean=mu_theta_3_posterior,sd=sqrt(sigma_theta_3_2_posterior))
    gamma_3_posterior=exp(rnorm(n=S,mean=mu_gamma_3_posterior,sd=sqrt(sigma_gamma_3_2_posterior)))
    Tau_3_posterior=rtruncnorm(n=S,a=(d_draws-Tau_star_3),b=d_draws,mean=d_draws-mu_tau_3_posterior,sd=sqrt(sigma_tau_3_2_posterior))
    mean_Y3=rep(theta_3_posterior,j) + rep(I_3_posterior,j)*rep(gamma_3_posterior,j)*(temp_t-Tau_3_posterior)*as.numeric(temp_t>Tau_3_posterior)
    
    prob_Y3_HCC=dnorm(temp_Y3,mean=mean_Y3,sd=sqrt(sigma_3_posterior))
    prob_Y3_HCC_1=array(prob_Y3_HCC,c(S,j))
    prob_Y3_HCC_2=apply(prob_Y3_HCC_1,1,prod,na.rm=TRUE)
    p_Y3_HCC[l]=mean(prob_Y3_HCC_2)
      
    l=l+1
  }
  #print(c(i,j,l))
}

#Get results
p_HCC=mean(data$D[(data$obs_number==1)]) #prior probability of HCC in training data
p_noHCC=1-p_HCC
data_validation["posterior_risk_3marker"]=(p_Y1_HCC*p_Y2_HCC*p_Y3_HCC*p_HCC)/(p_Y1_noHCC*p_Y2_noHCC*p_Y3_noHCC*p_noHCC)

########################################################################################################################################################################################################################
#Estimate ROC(0.1)
########################################################################################################################################################################################################################
specificity_fixed=0.9

cut_off_3marker=quantile(data_validation$posterior_risk_3marker[(data_validation$D==0)],probs=specificity_fixed,type=1)
at_least_one_positive=rep(NA,N-n_0)
for(i in 1:(N-n_0))
{
  subject_data=subset(data_validation,data_validation$ID==subject_ID[i+n_0]) 
  at_least_one_positive[i]=as.numeric(sum(subject_data$posterior_risk_3marker>cut_off_3marker,na.rm=T)>0)
}

round(mean(at_least_one_positive,)*100,2)
#77.55


