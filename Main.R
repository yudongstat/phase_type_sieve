rm(list = ls())
setwd("C:/Users/ydwang/Dropbox/Biometrika_revision/Code")
load("Covariate_mat.rdata")
source("Core_functions.R")
library(expm)
library(mapfit)
library(MCMCpack)
library(FAdist)

###################################################################################
######        Block 1. Data generation              #######
###################################################################################
sample_size<-100;# sample size
# generate reporting delay
lags=rweibull(sample_size,shape = 2,scale = 0.8)
# generate baseline lifetime
ts=rgamma(sample_size, shape = 2.5,scale = 0.9)
# generate baseline dropout time
tau<-rllog(sample_size,shape = 1/2.5,scale = log(3.5))
# generate censoring time
censor<-runif(n=sample_size,min = 0,max = 8)
# generate covariates
covariate_dim<-5 # dimension of covariate
# we generate the covariates by resampling from the covariates of the real data in Section 6 of our paper
covariate_all<-covariate
covariate_all<-covariate_all[sample(1:1252,1252),]
covariate<-covariate_all[sample(1:1252,size = sample_size,replace = TRUE),]

# set the true value for regression coefficients in the AFT model
true_beta<-c(0.24,  0.01,  0.42,  1, -0.25)
true_alpha<-c(0.18,  0.01,  0.15,  0.15, -0.06)

# generate the lifetime and dropout time
ts<-c(ts*exp(covariate%*%true_beta))
tau<-c(tau*exp(covariate%*%true_alpha))
original_data_list <- list(covariate, ts, tau, lags, censor)

# only those with min(T,D) + R < C will be observed
ids=which((pmin(ts,tau)+lags)<censor)
# reporting delays for reported patients
lags=lags[ids]
# observed lifetime/dropout time for reported patients
observed<-pmin(ts,tau)[ids]
# observed censoring indicator for reported patients
delta<-as.numeric(pmin(ts,tau)[ids]==ts[ids])
# covariates for reported patients
covariate_reported<-covariate[ids,]
# covariates for unreported patients
covariate_unreported<-covariate[-ids,]
# number of unreported patients
ncensor=sample_size-length(lags)
# the admin censoring time for unreported patients
censor<-censor[-ids]


# set the order of the phase-type distributions
m1<-ceiling((sample_size)^(1/4))
m2<-ceiling((sample_size)^(1/4))
m3<-ceiling((sample_size)^(1/4))



###################################################################################
######        Block 2. Tuning parameter selection              #######
###################################################################################
# this indicator controls whether to implement CV or not
# if TRUE, cross-validation will be implemented to determine (s1, s2, s3), which correspond to (\lambda, \overline{\lambda}, \tilde{\lambda}) in our paper
# if FALSE, we will randomly choose those tuning parameters
# be mindful that the CV procedure is time-consuming
CV_flag <- FALSE

# set candidate values for (s1, s2, s3), which correspond to (\lambda, \overline{\lambda}, \tilde{\lambda}) in our paper
s1_candidate <- seq(0.3,0.35,0.01)
s2_candidate <- seq(0.2,0.25,0.01)
s3_candidate <- seq(1.1,1.15,0.01)

if(CV_flag){
  # set the number of folds that you will use for CV
  fold <- 2
  # run CV; this is time-consuming!
  selected_s <- cross_validation(s1_candidate, s2_candidate, s3_candidate, fold, original_data_list)
  s1 <- selected_s[1]
  s2 <- selected_s[2]
  s3 <- selected_s[3]
}else{
  s1 <- sample(x = s1_candidate, size = 1)
  s2 <- sample(x = s2_candidate, size = 1)
  s3 <- sample(x = s3_candidate, size = 1)
}



###################################################################################
######        Block 3. Point estimation via EM algorithm              #######
###################################################################################

# define (p0ps, p1ps, p2ps), which corresponds to \Lambda, \overline{\Lambda}, \tilde{\Lambda} in our paper
# also define (nups, nu1ps, nu2ps), which corresponds to \xi, \overline{\xi}, \tilde{\xi} in our paper
initial_res <- initialization(s1, s2, s3, m1, m2, m3)
p0ps <- initial_res[[1]]
p1ps <- initial_res[[2]]
p2ps <- initial_res[[3]]
nups <- initial_res[[4]] 
nu1ps <- initial_res[[5]] 
nu2ps <- initial_res[[6]] 

# set initial values for old_pips, old_pi1ps, old_pi2ps, which corresponds to \pi, \overline{\pi}, and \tilde{\pi} in our paper
initial_pips<-rep(1/m1,m1)
initial_pi1ps<-rep(1/m2,m2)
initial_pi2ps<-rep(1/m3,m3)

# set initial values for beta and alpha
initial_beta<-rep(0,covariate_dim)
initial_alpha<-rep(0,covariate_dim)

# run the EM algorithm
EM_res <- EM_algo(observed, delta, lags, censor, initial_pips, initial_pi1ps, initial_pi2ps, initial_beta, initial_alpha, stopping_criteria = 1e-2, max_iter = 1000)

# point estimator
MLE_beta <- EM_res[[1]]
MLE_alpha <- EM_res[[2]]
pi_MLE1 <- EM_res[[3]]
pi_MLE2 <- EM_res[[4]]
pi_MLE3 <- EM_res[[5]]



###################################################################################
######        Block 4. Interval estimation              #######
###################################################################################
CI_res <- confidence_interval(original_data_list, p0ps, p1ps, p2ps, nups, nu1ps, nu2ps, EM_res)
CI_lower <- CI_res[[1]]
CI_upper <- CI_res[[2]]
