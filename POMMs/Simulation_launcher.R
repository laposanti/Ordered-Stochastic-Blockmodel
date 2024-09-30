

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

setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/")

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


folder_url <- "https://drive.google.com/drive/u/1/folders/14WD-QRXQbCvbPjqpHHBenrQuxwDXTzHd"
folder <- drive_get(as_id(folder_url))




################################################################################
#                        Choose the directory
################################################################################

#chosing where to save the files



print('Simulation study for fixed K, for K=3,4,5,6')

is.simulation=T


#data.directory

data_directory = "./Data/Sim_2_data//"

for(true_model in  c('SST','Simple')){
  true_model = 'SST'
  filenames <- list.files(pattern = true_model,path =data_directory)
  print(filenames) #data to be estimated
  

  
  choose_model_to_estimate = c('SST', 'WST', 'Simple')
  #-----------------------------------------------------------------------------
  # read the files in the selected folder, estimate the SST, the WST and the Simple model
  #-----------------------------------------------------------------------------
  
  for(file in 1:length(filenames)){
    file=1
    data_to_be_estimated = readRDS(paste0(data_directory,filenames[file]))
    data_to_be_estimated$data_plot
    data_to_be_estimated$sparsity
    stopifnot(data_to_be_estimated$model == true_model)
    recovery_capability = data_to_be_estimated$recovery_capability
    N_ij = data_to_be_estimated$N_ij
    n = nrow(N_ij)
    Y_ij = data_to_be_estimated$Y_ij
    ground_truth =data_to_be_estimated$ground_truth

    
    K= nrow(data_to_be_estimated$ground_truth$theta)
    data_description = paste0(true_model,K)
    print(paste0("True data--->", filenames[file]))
    
    ##############################################################################
    # Estimation: set the parameters of the estimation
    ##############################################################################
    
    n_chains = 4
    optimal_acceptance_rate_theta =.44
    optimal_acceptance_rate_mu = .234

    N_iter <- 40000  #number of iterations
    burnin <- 10000  #number of discarded iterations
    thin=5
    diag0.5=T
    seed=23
    # for(K in 3:10){
    K = data_to_be_estimated$ground_truth$K
    K_est = rep(K, n_chains) #number of clusters to fit
    
    is.simulation=T
    
    print(paste0("True data--->", filenames[file], "\n"))
    #-----------------------------------------------------------------------------
    # SST MODEL
    #-----------------------------------------------------------------------------
    if('SST' %in% choose_model_to_estimate){
      print(paste0("Estimation of the SST model, K=",K))
      print(paste0("Begin cycle at:", date()))
      
      
      
      est_model = 'SST'
      
      #setting up the chain hyperparameter
      
      #where to save the data
      saving_directory = paste0("./Results/MCMC_output/Fixed_K/Simulation_Nials'data//True_Model",true_model,"/")
      
      
      #Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
      power_posterior_apprach = F
      custom_init <- NA
      
      print(paste0("Estimation of the SST model, K=", K_est))
      print(paste0("Begin cycle at:", date(), "\n"))
      source("./model_auxiliary_functions/Functions_priorSST.R")
      source("./Metropolis_within_Gibbs_code_powerposterior.R")
      
      estimation_control <- list(z = 1,  theta = 1)
      
      chains_SST <- adaptive_MCMC_orderstats_powerposterior(Y_ij = Y_ij, N_ij = N_ij, 
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
                                                            thin=thin,diag0.5 = T)
      
      
      
      
      
      my_names <- paste0("chain", 1:n_chains)
      names(chains_SST)<- my_names 
      chains_SST[['recovery_level']] = recovery_capability
      
      
      my_filename = paste0('./Results/MCMC_output/Fixed_K/Simulation/Data_from',data_description, "_est_model",
                           est_model,"_Kest",paste(unlist(K_est),collapse = "_"),
                           'recovery_level',recovery_capability,'.rds')
      saveRDS(object = chains_SST, file = my_filename) 
      drive_put(media = my_filename, path = folder)
      
    }
    
    #-----------------------------------------------------------------------------
    # WST MODEL
    #-----------------------------------------------------------------------------
    
    if('WST' %in% choose_model_to_estimate){
      print(paste0("Estimation of the WST model, K=",K))
      print(paste0("Begin cycle at:",date()))
      #initializing each chain
      
      
      est_model = 'WST'
      
      #setting up the chain hyperparameter
      
      #where to save the data
      saving_directory = "./Results/"
      
      
      #Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
      power_posterior_apprach = F
      custom_init <- NA
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
                                                            thin = thin,diag0.5 = T)
      
      
      
      
      my_names <- paste0("chain", 1:n_chains)
      names(chains_WST)<-my_names 
      chains_WST[['recovery_level']] = recovery_capability
      my_filename = paste0('./Results/MCMC_output/Fixed_K/Simulation/Data_from',
                           data_description, "_est_model",
                           est_model,"_Kest",paste(unlist(K_est),collapse = "_"),
                           'recovery_level',recovery_capability,'.rds')
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
      
      #setting up the chain hyperparameter
      
      #where to save the data
      saving_directory = "./Results/"
      
      
      #Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
      power_posterior_apprach = F
      custom_init <- NA
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
                                                              thin=thin,diag0.5= F)
      my_names <- paste0("chain", 1:n_chains)
      names(chains_Simple)<- my_names 
      chains_Simple[['recovery_level']] = recovery_capability
      my_filename = paste0('./Results/MCMC_output/Fixed_K/Simulation/Data_from',
                           data_description, "_est_model",
                           est_model,"_Kest",paste(unlist(K_est),collapse = "_"),
                           'recovery_level',recovery_capability,'.rds')
      saveRDS(object = chains_Simple, file = my_filename) 
      
      drive_put(media = my_filename, path = folder)
    }
  }

}




