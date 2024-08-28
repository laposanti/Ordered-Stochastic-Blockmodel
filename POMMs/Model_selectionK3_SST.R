
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
library(doRNG)
library(googledrive)
#setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/")
source("./model_auxiliary_functions/Functions_priorSST.R")
source("./Metropolis_within_Gibbs_code.R")
source("./Metropolis_within_Gibbs_code_powerposterior.R")




################################################################################
#                        step 1: upload the data
################################################################################


#where the data are stored
data_wd<- "./Data/Sim1_data/"
data_description = 'SST3'
filenames <- list.files(pattern = paste0(data_description),path = data_wd)
print(filenames)
data_to_be_estimated <- readRDS(paste0(data_wd, "/", filenames[1]))
recovery_capability = data_to_be_estimated$recovery_capability
N_ij <- data_to_be_estimated$N_ij
n <- nrow(N_ij)
Y_ij <- data_to_be_estimated$Y_ij
K <- data_to_be_estimated$ground_truth$K
ground_truth <- data_to_be_estimated$ground_truth


subject = "lapo.santi@ucdconnect.ie"
service_account_key = "./sonic-426715-75af23aca274.json"
googledrive::drive_deauth()
googledrive::drive_auth_configure(path = "./client_secret_573831164304-jqqj3i5mhvubbkkuifvtgkfsut8lse3g.apps.googleusercontent.com.json")
googledrive::drive_auth(email = subject)
# 
# 
# 2
# # Get the folder (if you already have it) or specify the path where you want to upload
folder_url <- "https://drive.google.com/drive/u/1/folders/1QLqA5DfE1LSfqq7GJqDy3u5O2K7vwc8N"
folder <- drive_get(as_id(folder_url))




################################################################################
# Decide for how many Ks we want to compute the marginal posterior
################################################################################

print(paste0("True data--->", filenames[1], "\n"))
is.simulation=T
choose_model_to_estimate = c('SST',"WST","Simple")

#setting up the chain hyperparameter
N_iter <- 60000  #number of iterations
burnin <- 30000 #number of discarded iterations

K_est = list(2,3,4,5,6,7,8,9,10) #number of clusters to fit
#where to save the data
saving_directory = "./Results/MCMC_output/model_choice/WAIC_method/K3_true//"

# Check if the directory exists
if (!dir.exists(saving_directory)) {
  # If the directory doesn't exist, create it
  dir.create(saving_directory, recursive = TRUE)
  message("Directory created.")
} else {
  message("Directory already exists.")
}
#Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
power_posterior_apprach = F
custom_init <- NA


if('SST' %in% choose_model_to_estimate){
  
  

  


  print(paste0("Estimation of the SST model, K=", K_est))
  print(paste0("Begin cycle at:", date(), "\n"))
  
  est_model = 'SST'
  seed=23
  
  estimation_control <- list(z = 1, sigma_squared = 0, mu_vec = 1, K = 0, theta = 0)
  
  chains<- adaptive_MCMC_orderstats_powerposterior(Y_ij = Y_ij, N_ij = N_ij, 
                                                   saving_directory = saving_directory,
                                                   estimation_control = estimation_control,
                                                   burnin = burnin,
                                                   ground_truth = ground_truth,
                                                   n = n, N_iter = N_iter, 
                                                   K_est = K_est,data_description = data_description,
                                                   seed = seed, 
                                                   model = est_model, 
                                                   custom_init = custom_init,
                                                   power_posterior_apprach = power_posterior_apprach,
                                                   thin=thin)
  names(chains) = paste0('chain',unlist(K_est))
  
  
  chains[['recovery_level']] = recovery_capability
  
  my_filename = paste0(saving_directory,data_description, "_est_model",
                       est_model,"_Kest",paste(unlist(K_est),collapse = "_"),
                       'recovery_level',recovery_capability,'.rds')
  saveRDS(object = chains, file = my_filename) 
  drive_put(media = my_filename, path = folder)
  
}

if('WST' %in% choose_model_to_estimate){
  print(paste0("Estimation of the WST model, K=",K))
  print(paste0("Begin cycle at:",date()))
  #initializing each chain
  
  
  est_model = 'WST'
  
  #Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
  power_posterior_apprach = F
  custom_init <- NA
  print(paste0("Estimation of the WST model, K=", K_est))
  print(paste0("Begin cycle at:", date(), "\n"))
  estimation_control <- list(z = 1, sigma_squared = 0, mu_vec = 0 ,K = 0, theta = 1)
  
  chains_WST <- adaptive_MCMC_orderstats_powerposterior(Y_ij = Y_ij, N_ij = N_ij,
                                                        saving_directory = saving_directory,
                                                        estimation_control = estimation_control,
                                                        burnin = burnin,
                                                        ground_truth = ground_truth,
                                                        n = n, N_iter = N_iter, 
                                                        K_est = K_est,data_description = data_description,
                                                        seed = seed, 
                                                        model = est_model, 
                                                        custom_init = custom_init,
                                                        power_posterior_apprach = power_posterior_apprach,thin = thin)
  
  
  
  
  names(chains_WST) = paste0('chain',unlist(K_est))
  chains_WST[['recovery_level']] = recovery_capability
  
  my_filename = paste0(saving_directory,
                       data_description, "_est_model",
                       est_model,"_Kest",paste(unlist(K_est),collapse = "_"),
                       'recovery_level',
                       recovery_capability,'.rds')
  saveRDS(object = chains_WST, file = my_filename) 
  
  drive_put(media = my_filename, path = folder)
  
}

#-----------------------------------------------------------------------------
# Simple model
#-----------------------------------------------------------------------------


if('Simple' %in% choose_model_to_estimate){
  
  print(paste0("Estimation of Simple model, K=",K))
  print(paste0("Begin cycle at:",date()))
  
  
  est_model = 'Simple'
  
  
  #Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
  power_posterior_apprach = F
  custom_init <- NA
  estimation_control = list(z = 1,sigma_squared=0, mu_vec=0,K=0,theta=1)
  
  chains_Simple = adaptive_MCMC_orderstats_powerposterior(Y_ij = Y_ij, N_ij = N_ij,
                                                          saving_directory = saving_directory,
                                                          estimation_control = estimation_control,
                                                          burnin = burnin,
                                                          ground_truth = ground_truth,
                                                          n = n, N_iter = N_iter, 
                                                          K_est = K_est,data_description = data_description,
                                                          seed = seed, 
                                                          model = est_model, 
                                                          custom_init = custom_init,
                                                          power_posterior_apprach = power_posterior_apprach,
                                                          thin=thin)
  names(chains_Simple) = paste0('chain',unlist(K_est))
  chains_Simple[['recovery_level']] = recovery_capability
  my_filename = paste0(saving_directory,
                       data_description, "_est_model",
                       est_model,"_Kest",paste(unlist(K_est),collapse = "_"),
                       'recovery_level',recovery_capability,'.rds')
  saveRDS(object = chains_Simple, file = my_filename) 
  
  drive_put(media = my_filename, path = folder)
}