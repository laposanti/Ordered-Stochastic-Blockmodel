
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
source("./Metropolis_within_Gibbs_code_powerposterior.R")


subject = "lapo.santi@ucdconnect.ie"
service_account_key = "./sonic-426715-75af23aca274.json"
googledrive::drive_deauth()
googledrive::drive_auth_configure(path = "./client_secret_573831164304-jqqj3i5mhvubbkkuifvtgkfsut8lse3g.apps.googleusercontent.com.json")
googledrive::drive_auth(email = subject)
# 
# 
# 2
# # Get the folder (if you already have it) or specify the path where you want to upload






################################################################################
#                        step 1: upload the data
################################################################################


#where the data are stored
data_wd<- "./Data/Sim1_data/"
data_description = 'SST5'
filenames <- list.files(pattern = paste0(data_description),path = data_wd)
data_to_be_estimated <- readRDS(paste0(data_wd, "/", filenames[1]))
recovery_capability = data_to_be_estimated$recovery_capability
N_ij <- data_to_be_estimated$N_ij
n <- nrow(N_ij)
Y_ij <- data_to_be_estimated$Y_ij
K <- data_to_be_estimated$ground_truth$K
ground_truth <- data_to_be_estimated$ground_truth

choose_model_to_estimate = c('SST', 'WST','Simple')

################################################################################
# Decide for how many Ks we want to compute the marginal posterior
################################################################################

print(paste0("True data--->", filenames[1], "\n"))
is.simulation=T




optimal_acceptance_rate_theta =.44
optimal_acceptance_rate_mu = .234
seed = 23
N_iter <- 70000 #number of iterations
burnin <- 20000 #number of discarded iterations
thin = 5


K_est = list(2,3,4,5,6,7,8,9,10) #number of clusters to fit

#Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
power_posterior_apprach = F
custom_init <- NA

#where to save the data
saving_directory = './Results/MCMC_output/model_choice/WAIC_method/K5_true/'

# Check if the directory exists
if (!dir.exists(saving_directory)) {
  # If the directory doesn't exist, create it
  dir.create(saving_directory, recursive = TRUE)
  message("Directory created.")
} else {
  message("Directory already exists.")
}

power_posterior_apprach = F
custom_init = NA
for(diag0.5 in c(T,F)){
  if(diag0.5==T){
    #main diagonal fixed to 0.5
    folder_url <- "https://drive.google.com/drive/u/1/folders/1cneKgDKZ8ZdxNZ4_9-QACyoYDgl9c_9K"
  }else{
    folder_url <- "https://drive.google.com/drive/u/1/folders/1WMb0GZuW1Je4crjh3Fv5TtjwhQF-k24G"
    
  }
  folder <- drive_get(as_id(folder_url))
  if('SST' %in% choose_model_to_estimate){
    
    
    
    
    print(paste0("Estimation of the SST model, K=", K_est))
    print(paste0("Begin cycle at:", date(), "\n"))
    est_model = 'SST'
    
    seed=23
    
    estimation_control <- list(z = 1, theta = 1)
    
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
                                                      power_posterior_apprach = power_posterior_apprach,
                                                      thin =thin,
                                                      diag0.5 = diag0.5)
    names(chains) = paste0('chain',unlist(K_est))
    chains[['recovery_level']] = recovery_capability
    
    my_filename = paste0(saving_directory, 'Data_from',data_description, "_est_model",
                         est_model,"_Kest",paste(unlist(K_est),collapse = "_"),
                         'recovery_level',
                         recovery_capability,'.rds')
    saveRDS(object = chains, file = my_filename) 
    drive_put(media = my_filename, path = folder)
    
  }
  
  if('WST' %in% choose_model_to_estimate){
    print(paste0("Estimation of the WST model, K=",K))
    print(paste0("Begin cycle at:",date()))
    #initializing each chain
    
    
    est_model = 'WST'
    
    #Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
    
    print(paste0("Estimation of the WST model, K=", K_est))
    print(paste0("Begin cycle at:", date(), "\n"))
    estimation_control <- list(z = 1, theta = 1)
    
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
                                                          power_posterior_apprach = power_posterior_apprach,
                                                          thin = thin,
                                                          diag0.5 = diag0.5)
    
    
    
    
    names(chains_WST) = paste0('chain',unlist(K_est))
    chains_WST[['recovery_level']] = recovery_capability
    
    my_filename = paste0(saving_directory,"Data_from",
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
    
    estimation_control = list(z = 1,theta=1)
    
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
                                                            thin=thin,diag0.5 = F)
    names(chains_Simple) = paste0('chain',unlist(K_est))
    chains_Simple[['recovery_level']] = recovery_capability
    my_filename = paste0(saving_directory,"Data_from",
                         data_description, "_est_model",
                         est_model,"_Kest",paste(unlist(K_est),collapse = "_"),
                         'recovery_level',recovery_capability,'.rds')
    saveRDS(object = chains_Simple, file = my_filename) 
    
    drive_put(media = my_filename, path = folder)
  }
}
