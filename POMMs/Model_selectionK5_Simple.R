
library(mcclust)
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
library(ggside)
library(parallel)
library(truncnorm)
library(label.switching)
library(doRNG)

#setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/")
source("./model_auxiliary_functions/Functions_priorSST.R")
source("./Metropolis_within_Gibbs_code.R")
source("./model_auxiliary_functions/MCMC_functions.R")
source("./Metropolis_within_Gibbs_code_powerposterior.R")




################################################################################
#                        step 1: upload the data
################################################################################


#where the data are stored
data_wd<- "./Data/power_posterior_data/"
data_description = 'SST_K5'
filenames <- list.files(pattern = paste0(data_description),path = data_wd)
data_to_be_estimated <- readRDS(paste0(data_wd, "/", filenames[3]))
N_ij <- data_to_be_estimated$N_ij
n <- nrow(N_ij)
Y_ij <- data_to_be_estimated$Y_ij
K <- data_to_be_estimated$ground_truth$K
ground_truth <- data_to_be_estimated$ground_truth



################################################################################
# Decide for how many Ks we want to compute the marginal posterior
################################################################################

print(paste0("True data--->", filenames[3], "\n"))
is.simulation=T

est_model = 'Simple'

#setting up the chain hyperparameter
N_iter <- 120000 #number of iterations
burnin <- 80000  #number of discarded iterations

K_est = list(2,3,4,5,6,7,8) #number of clusters to fit


#where to save the data
saving_directory = "./Results/"


#Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
power_posterior_apprach = T
custom_init <- NA

print(paste0("Estimation of the", est_model," model, K=", K_est))
print(paste0("Begin cycle at:", date(), "\n"))


seed=23

estimation_control <- list(z = 1, sigma_squared = 0, mu_vec = 0, K = 0, theta = 1)

chains <- adaptive_MCMC_orderstats_powerposterior(Y_ij = Y_ij, N_ij = N_ij,
                                                  saving_directory = saving_directory,
                                                  estimation_control = estimation_control,
                                                  burnin = burnin,
                                                  ground_truth = ground_truth,
                                                  n = n, N_iter = N_iter, 
                                                  K_est = K_est,data_description = data_description,
                                                  seed = seed, 
                                                  model = est_model, 
                                                  custom_init = custom_init,
                                                  power_posterior_apprach = T)
names(chains) = paste0('chain',unlist(K_est))





