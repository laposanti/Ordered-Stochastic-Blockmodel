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

# a short description of the data
data_description = 'Citations_application'

Y_ij <- read.csv("./Data/Citations_application/cross-citation-matrix.csv",header = T,row.names = 1)
diag(Y_ij) = 0
N_ij= matrix(0,47,47) +Y_ij*upper.tri(Y_ij)+t(Y_ij)*upper.tri(Y_ij)+Y_ij*lower.tri(Y_ij)+t(Y_ij)*lower.tri(Y_ij)

Y_ij = as.matrix(Y_ij)
N_ij = as.matrix(N_ij)

n <- nrow(N_ij)


################################################################################
# Decide for how many Ks we want to compute the marginal posterior
################################################################################

#--------------------------------
#which model to fit
est_model = 'Simple'

#which parameters to estimate: 1 estimated, 0 fixed. if 0 is selected, provide some initial quantity
estimation_control <- list(z = 1, sigma_squared = 0, mu_vec = 0, K = 0, theta = 1)
custom_init <- NA #whether you want to initialize the parameters manually
ground_truth = NA

#----------------------------------
#setting chains' features

N_iter <- 50000 #number of iterations
burnin <- 20000  #number of discarded iterations #number of iterations
K_est = list(2,3,4,5,6,7,8) #K values to be fitted
power_posterior_apprach = T #Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
seed=23 #the seed for the chains
#----------------------------------
#where to save the data
saving_directory = "./Results/"


print(paste0("Estimation of the WST model, K=", K_est))
print(paste0("Begin cycle at:", date(), "\n"))


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




