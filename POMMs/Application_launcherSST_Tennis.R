

library(doFuture)
library(progressr)
library(beepr)
library(foreach)
library(doParallel)
library(tidyverse)
library(EnvStats)
library(truncnorm)
library(dplyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(truncnorm)
library(doRNG)




source("./model_auxiliary_functions/Functions_priorSST.R")
source("./Metropolis_within_Gibbs_code.R")
source("./model_auxiliary_functions/MCMC_functions.R")
source("./Metropolis_within_Gibbs_code_powerposterior.R")



###############################################################################
# uploading data
###############################################################################

data_description = 'Tennis_application'

Y_ij <- read.table("./Data/Tennis application/Y_ij.csv",header  = F,row.names = 1,sep = ",")
N_ij <- read.table("./Data/Tennis application/N_ij.csv",header  = F,row.names = 1,sep = ",")

Y_ij = as.matrix(Y_ij)
N_ij = as.matrix(N_ij)

n <- nrow(N_ij)



################################################################################
# Decide for how many Ks we want to compute the marginal posterior
################################################################################



est_model = 'SST'

#setting up the chain hyperparameter
N_iter <- 120000 #number of iterations
burnin <- 80000  #number of discarded iterations

K_est = list(2,3,4,5,6,7,8) #K values to be fitted

#where to save the data
saving_directory = "./Results/"


#Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
power_posterior_apprach = T 

#whether you want to initialize the parameters manually
custom_init <- NA
ground_truth = NA
seed=23
estimation_control <- list(z = 1, sigma_squared = 0, mu_vec = 1, K = 0, theta = 1)


print(paste0("Estimation of the ", est_model," model, K=", K_est))
print(paste0("Begin cycle at:", date(), "\n"))



################################################################################
# Decide for how many Ks we want to compute the marginal posterior
################################################################################

chains <- adaptive_MCMC_orderstats_powerposterior(Y_ij = Y_ij, N_ij = N_ij,
                                                  saving_directory = saving_directory,
                                                  estimation_control = estimation_control,
                                                  burnin = burnin,
                                                  ground_truth = ground_truth,
                                                  n = n, N_iter = N_iter, 
                                                  K_est = K_est,
                                                  data_description = data_description,
                                                  seed = seed, 
                                                  model = est_model, 
                                                  custom_init = custom_init,
                                                  power_posterior_apprach = power_posterior_apprach)
names(chains) = paste0('chain',unlist(K_est))
