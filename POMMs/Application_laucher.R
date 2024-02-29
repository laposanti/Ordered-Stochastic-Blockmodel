



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

source("/Users/lapo_santi/Desktop/Nial/oldmaterial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Metropolis_within_Gibbs_code.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/model_auxiliary_functions/MCMC_functions.R")





################################################################################
#                         Set up parameters for the simulation study
################################################################################






true_model = 'Tennis_data'

###############################################################################
# uploading data
###############################################################################



#chosing where to save the files
setwd("/Users/lapo_santi/Desktop/Nial/MCMC_results/application_13Feb2024//raw/")


 #data to be estimated


K_values <- c(3,4,5,6)  # Range of K values to explore


choose_model_to_estimate = c('SST', 'WST','Simple')
#-----------------------------------------------------------------------------
# read the files in the selected folder, estimate the SST, the WST and the Simple model
#-----------------------------------------------------------------------------

for(k_th in K_values){
  

  
  Y_ij <- read.table("/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/rawdata/Y_ij.csv",header  = F,row.names = 1,sep = ",")
  N_ij <- read.table("/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/rawdata/N_ij.csv",header  = F,row.names = 1,sep = ",")
  
  Y_ij = as.matrix(Y_ij)
  N_ij = as.matrix(N_ij)
  
  n = nrow(N_ij)
  K = list(k_th,k_th,k_th,k_th)
  ground_truth = NA
  
    
  ##############################################################################
  # Estimation: set the parameters of the estimation
  ##############################################################################
  
  n_chains = 4
  optimal_acceptance_rate =.22
  N_iter= 110000
  chains_seeds = list(20,09,97,2024)
  
  #-----------------------------------------------------------------------------
  # SST MODEL
  #-----------------------------------------------------------------------------
  if('SST' %in% choose_model_to_estimate){
    print(paste0("Estimation of the SST model, K=",K))
    print(paste0("Begin cycle at:",date()))
    #initializing each chain
    
    
    
    estimation_control = list(z = 1,sigma_squared=0, mu_vec=1,K=0,P=1)
    
    
    
    chains_SST = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij , 
                                          estimation_control = estimation_control, 
                                          ground_truth = ground_truth, 
                                          n = n, N_iter = N_iter,n_chains = n_chains, 
                                          optimal_acceptance_rate=optimal_acceptance_rate, K = K,
                                          seed = chains_seeds, model = 'SST')
    
    
    my_names <- paste0("chain", 1:n_chains)
    names(chains_SST)<-my_names 
    
    
    
    filename_SST <- paste0("True_Model",true_model,"Est_model_SST","_N", n,"_K", K[[1]],".RDS")
    saveRDS(chains_SST, file = filename_SST) #saving results
    beep("coin")
  }
  
  #-----------------------------------------------------------------------------
  # WST MODEL
  #-----------------------------------------------------------------------------
  
  if('WST' %in% choose_model_to_estimate){
    print(paste0("Estimation of the WST model, K=",K))
    print(paste0("Begin cycle at:",date()))
    #initializing each chain
    
    estimation_control = list(z = 1,sigma_squared=1, mu_vec=1,K=0,P=1)
    
    
    chains_WST = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij , 
                                          estimation_control = estimation_control, 
                                          ground_truth = ground_truth, 
                                          n = n, N_iter = N_iter,n_chains = n_chains, 
                                          optimal_acceptance_rate=optimal_acceptance_rate, K = K,
                                          seed = chains_seeds, model = 'WST')
    my_names <- paste0("chain", 1:n_chains)
    names(chains_WST)<-my_names 
    
    filename_WST <- paste0("True_Model",true_model,"Est_model_WST","_N", n,"_K", K[[1]],".RDS")
    saveRDS(chains_WST, file = filename_WST) #saving results
    beep("coin")
    
  }
  #-----------------------------------------------------------------------------
  # Simple model
  #-----------------------------------------------------------------------------
  
  
  if('Simple' %in% choose_model_to_estimate){
    print(paste0("Estimation of Simple model, K=",K))
    print(paste0("Begin cycle at:",date()))
    
    
    estimation_control = list(z = 1,sigma_squared=0, mu_vec=0,K=0,P=1)
    
    
    
    chains_Simple = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij , 
                                             estimation_control = estimation_control, 
                                             ground_truth = ground_truth, 
                                             n = n, N_iter = N_iter,n_chains = n_chains, 
                                             optimal_acceptance_rate=optimal_acceptance_rate, K = K,
                                             seed = chains_seeds, model = 'Simple')
    my_names <- paste0("chain", 1:n_chains)
    names(chains_Simple)<- my_names 
    filename_Simple <- paste0("True_Model",true_model,"Est_model_Simple","_N", n,"_K", K[[1]],".RDS")
    saveRDS(chains_Simple, file = filename_Simple) #saving results
    beep("coin")
  }
}


