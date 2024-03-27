

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




################################################################################
#                         Set up parameters for the simulation study
################################################################################





#choose between citation exchange data and tennis data
#citation data::: set true_model =  "Citation_data"
#tennis data::: tennis data = 'Tennis_data'
for(application in c('Tennis_data', "Citation_data")){
  true_model = application
  
  ###############################################################################
  # uploading data
  ###############################################################################
  
  if(true_model == 'Tennis_data'){
    Y_ij <- read.table("./Data/Tennis application/Y_ij.csv",header  = F,row.names = 1,sep = ",")
    N_ij <- read.table("./Data/Tennis application/N_ij.csv",header  = F,row.names = 1,sep = ",")
    
    Y_ij = as.matrix(Y_ij)
    N_ij = as.matrix(N_ij)
  }else if(true_model == 'Citation_data'){
    Y_ij <- read.csv("./Data/Citations_application/cross-citation-matrix.csv",header = T,row.names = 1)
    diag(Y_ij) = 0
    N_ij= matrix(0,47,47) +Y_ij*upper.tri(Y_ij)+t(Y_ij)*upper.tri(Y_ij)+Y_ij*lower.tri(Y_ij)+t(Y_ij)*lower.tri(Y_ij)
    
    Y_ij = as.matrix(Y_ij)
    N_ij = as.matrix(N_ij)
  }
  
  
  #chosing where to save the files depending on which model you are estimating
  
  
  #data to be estimated
  
  
  choose_model_to_estimate = c('WST')
  #-----------------------------------------------------------------------------
  # read the files in the selected folder, estimate the SST, the WST and the Simple model
  #-----------------------------------------------------------------------------
  
  
  n = nrow(N_ij)
  
  ground_truth = NA
  
  
  ##############################################################################
  # Estimation: set the parameters of the estimation
  ##############################################################################
  optimal_acceptance_rate_P =.44
  optimal_acceptance_rate_mu = .234
  N_iter <- 120000
  n_temperatures = 50
  #-----------------------------------------------------------------------------
  # SST MODEL
  #-----------------------------------------------------------------------------
  
  if('WST' %in% choose_model_to_estimate){
    
    
    
    K_est = list(2,3,4,5,6,7)
    saving_directories = list()
    chains_seeds <- list()
    for(i in 1:length(K_est)){
      saving_directories[[i]]<- paste0("./results/application/",application,"/WST/K", K_est[[i]],"//")
      chains_seeds[[i]] = 20+i
    }
    
    
    custom_init <- NA
    
    n_chains <- 6
    true_model <- application
    
    
    cat(paste0("Estimation of the WST model, K=", K_est, "\n"))
    cat(paste0("Begin cycle at:", date(), "\n"))
    
    estimation_control <- list(z = 1, sigma_squared =1, mu_vec = 1, K = 0, P = 1)
    
    chains <- adaptive_MCMC_orderstats_powerposterior(Y_ij = Y_ij, N_ij = N_ij,n_temperatures = n_temperatures,
                                                      saving_directory = saving_directories,
                                                      estimation_control = estimation_control,
                                                      ground_truth = ground_truth,
                                                      n = n, N_iter = N_iter, n_chains = n_chains,
                                                      optimal_acceptance_rate_P = optimal_acceptance_rate_P, 
                                                      optimal_acceptance_rate_mu = optimal_acceptance_rate_mu,
                                                      K = K_est,true_model = true_model,
                                                      seed = chains_seeds, est_model = 'WST', custom_init = custom_init)
    names(chains) = paste0('chain',unlist(K_est))
    
    
    
    
  }
  
  
}

