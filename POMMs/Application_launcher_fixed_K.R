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
library(googledrive)


#setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/")

source("./model_auxiliary_functions/Functions_priorSST.R")
source("./Metropolis_within_Gibbs_code_powerposterior.R")
source("./model_auxiliary_functions/MCMC_functions.R")

# Define the path to your service account key file


service_account_key <- "./client_secret_573831164304-jqqj3i5mhvubbkkuifvtgkfsut8lse3g.apps.googleusercontent.com.json"
t <- gargle::credentials_service_account(
  path = service_account_key,
  subject = "lapo.santi@ucdconnect.ie"
)
googledrive::drive_auth(token = t)
# Specify the Google Drive folder URL or ID
2
# Get the folder (if you already have it) or specify the path where you want to upload
folder_url <- "https://drive.google.com/drive/u/1/folders/1V-lQDh0DCWSx57YJ1hHf7ebwd6UinY6Z"
folder <- drive_get(as_id(folder_url))
################################################################################
#                         Set up parameters for the simulation study
################################################################################





#choose between citation exchange data and tennis data
#citation data::: set true_model =  "Citation_data"
#tennis data::: tennis data = 'Tennis_data'

for(data_description in c("Tennis_data","Citation_data")){

  
  ###############################################################################
  # uploading data
  ###############################################################################
  
  if(data_description == 'Tennis_data'){
    Y_ij <- readRDS("./Data/Tennis application/Y_new.rds")
    N_ij <- readRDS("./Data/Tennis application/N_new.rds")
    
    Y_ij = as.matrix(Y_ij)
    N_ij = as.matrix(N_ij)
    
  }else if(data_description == 'Citation_data'){
    Y_ij <- read.csv("./Data/Citations_application/cross-citation-matrix.csv",header = T,row.names = 1)
    diag(Y_ij) = 0
    N_ij= matrix(0,47,47) +Y_ij*upper.tri(Y_ij)+t(Y_ij)*upper.tri(Y_ij)+Y_ij*lower.tri(Y_ij)+t(Y_ij)*lower.tri(Y_ij)
    
    Y_ij = as.matrix(Y_ij)
    N_ij = as.matrix(N_ij)
  }
  
  
  #chosing where to save the files depending on which model you are estimating
  
  K_values <- c(3:8)  # Range of K values to explore
  
  print(paste0('Fitting now:' , data_description))
  
  choose_model_to_estimate = c("SST")
  #-----------------------------------------------------------------------------
  # read the files in the selected folder, estimate the SST, the WST and the Simple model
  #-----------------------------------------------------------------------------
  
  for(k_th in K_values){

    
    ##############################################################################
    # Estimation: set the parameters of the estimation
    ##############################################################################
    ground_truth = NA
    n = nrow(N_ij)

    power_posterior_apprach=F
    n_chains = 4
    optimal_acceptance_rate_theta =.44
    optimal_acceptance_rate_mu = .234
    seed=20
    N_iter <- 800 #number of iterations
    burnin <- 200 #number of discarded iterations
    
    
    K_est = rep(k_th, n_chains) #number of clusters to fit
    #-----------------------------------------------------------------------------
    # SST MODEL
    #-----------------------------------------------------------------------------
    if('SST' %in% choose_model_to_estimate){
      print(paste0("Estimation of the SST model, K=",k_th))
      print(paste0("Begin cycle at:", date()))
      
      
      
      est_model = 'SST'
      
      #where to save the data
      saving_directory = "./Results/"
      
      
      #Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
      power_posterior_apprach = F
      custom_init <- NA
      
      print(paste0("Estimation of the SST model, K=", K_est))
      print(paste0("Begin cycle at:", date(), "\n"))
      
      
      estimation_control <- list(z = 1, mu_vec = 1 ,K = 0, theta = 0)
      
      
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
                                                            power_posterior_apprach = power_posterior_apprach,thin=15)
      
      my_names <- paste0("chain", 1:n_chains)
      names(chains_SST)<- my_names 
      
      my_filename = paste0('./Results/MCMC_output/Fixed_K/Application/Data_from',data_description, "_est_model",est_model,"_Kest",K_est[[1]],'.rds')
      saveRDS(object = chains_SST, file = my_filename) 

      drive_put(media = my_filename, path = folder)
      
      
      
      
    }
    
    #-----------------------------------------------------------------------------
    # WST MODEL
    #-----------------------------------------------------------------------------
    
    if('WST' %in% choose_model_to_estimate){
      print(paste0("Estimation of the WST model, K=",k_th))
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
                                                            power_posterior_apprach = power_posterior_apprach,thin=15)
      
      
      
      
      my_names <- paste0("chain", 1:n_chains)
      names(chains_WST)<-my_names 
      my_filename = paste0('./Results/MCMC_output/Fixed_K/Application/Data_from',data_description, "_est_model",est_model,"_Kest",K_est[[1]],'.rds')
      saveRDS(object = chains_WST, file = my_filename) 
      
      
      drive_put(media = my_filename, path = folder)
      
      
      
      
    }
    
    #-----------------------------------------------------------------------------
    # Simple model
    #-----------------------------------------------------------------------------
    
    
    if('Simple' %in% choose_model_to_estimate){
      
      print(paste0("Estimation of Simple model, K=",k_th))
      print(paste0("Begin cycle at:",date()))
      
      
      est_model = 'Simple'
      
      #setting up the chain hyperparameter
      
      #where to save the data
      saving_directory = "./Results/"
      
      
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
                                                              power_posterior_apprach = power_posterior_apprach,thin=15)
      my_names <- paste0("chain", 1:n_chains)
      names(chains_Simple)<- my_names 
      my_filename = paste0('./Results/MCMC_output/Fixed_K/Application/Data_from',data_description, "_est_model",est_model,"_Kest",K_est[[1]],'.rds')
      saveRDS(object = chains_Simple, file = my_filename) 
      
      drive_put(media = my_filename, path = folder)
      
      
    }
  
  
  }
  
  
}

