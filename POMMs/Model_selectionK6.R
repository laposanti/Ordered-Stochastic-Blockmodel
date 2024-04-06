
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

is.simulation=T

true_model = 'SST'
est_model = 'SST'
#where the data are stored
#where the data are stored
data_wd<- "./Data/power_posterior_data/"

filenames <- list.files(pattern = paste0(true_model),path = data_wd)

choose_model_to_estimate=c('SST','WST','Simple')

#choose which data you want to use
# Choose which data you want to use
file=4
data_to_be_estimated <- readRDS(paste0(data_wd, "/", filenames[file]))
N_ij <- data_to_be_estimated$N_ij
n <- nrow(N_ij)
Y_ij <- data_to_be_estimated$Y_ij
K <- data_to_be_estimated$ground_truth$K
ground_truth <- data_to_be_estimated$ground_truth

print(paste0("True data--->", filenames[file], "\n"))

################################################################################
# Decide for how many Ks we want to compute the marginal posterior
################################################################################
n_temperatures=50
optimal_acceptance_rate_theta =.44
optimal_acceptance_rate_mu = .234
N_iter <- 120000

if('SST' %in% choose_model_to_estimate){
  
  
  K_est = list(2,3,4,5,6,7)
  saving_directories = list()
  chains_seeds <- list()
  for(i in 1:length(K_est)){
    saving_directories[[i]]<-  paste0("./results/model_selection/K", file + 2, "true/raw/K", K_est[[i]],"//")
    chains_seeds[[i]] = 20+i
  }
  
  
  print(paste0("Estimation of the SST model, K=", K_est))
  print(paste0("Begin cycle at:", date(), "\n"))
  custom_init <- NA
  
  n_chains <- 6
  estimation_control <- list(z = 1, sigma_squared = 0, mu_vec = 1, K = 0, theta = 1)
  
  chains <- adaptive_MCMC_orderstats_powerposterior(Y_ij = Y_ij, N_ij = N_ij,n_temperatures = n_temperatures,saving_directory = saving_directories,
                                                    estimation_control = estimation_control,
                                                    ground_truth = ground_truth,
                                                    n = n, N_iter = N_iter, n_chains = n_chains, 
                                                    optimal_acceptance_rate_theta  = optimal_acceptance_rate_theta, 
                                                    optimal_acceptance_rate_mu = optimal_acceptance_rate_mu,
                                                    K = K_est,true_model = true_model,
                                                    seed = chains_seeds,model = 'SST', custom_init = custom_init)
  names(chains) = paste0('chain',unlist(K_est))
  
  
  
  
}
