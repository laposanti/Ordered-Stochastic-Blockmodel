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


# Define the path to your service account key file

subject = "lapo.santi@ucdconnect.ie"
service_account_key = "./sonic-426715-75af23aca274.json"
googledrive::drive_deauth()
googledrive::drive_auth_configure(path = "./client_secret_573831164304-jqqj3i5mhvubbkkuifvtgkfsut8lse3g.apps.googleusercontent.com.json")
googledrive::drive_auth(email = subject)
# 
# 
# 2
#folder for application_fixed_K

################################################################################
#                         Set up parameters for the simulation study
################################################################################





#choose between citation exchange data and tennis data
#citation data::: set true_model =  "Citation_data"
#tennis data::: tennis data = 'Tennis_data'

for(data_description in c("Citation_data")){
  
  ###############################################################################
  # uploading data
  ###############################################################################
  
  if(data_description == 'Tennis_data'){
    Y_ij <- readRDS("./Data/Tennis application/Y_new.rds")
    N_ij <- readRDS("./Data/Tennis application/N_new.rds")
    
    Y_ij = as.matrix(Y_ij)
    N_ij = as.matrix(N_ij)
    
    
    
    
    
    pl_df = data.frame(x = rownames(N_ij), y = rowSums(N_ij))
    pl_df%>%
      ggplot(aes(x = x, y=y))+
      geom_col()+
      theme(axis.text.x = element_text(angle=90))
    
    
    
    
    
    
    
  }else if(data_description == 'Citation_data'){
    Y_ij <- read.csv("./Data/Citations_application/cross-citation-matrix.csv",header = T,row.names = 1)
    diag(Y_ij) = 0
    N_ij= matrix(0,47,47) +Y_ij*upper.tri(Y_ij)+t(Y_ij)*upper.tri(Y_ij)+Y_ij*lower.tri(Y_ij)+t(Y_ij)*lower.tri(Y_ij)
    
    Y_ij = as.matrix(Y_ij)
    N_ij = as.matrix(N_ij)
    
    
    
  }
  

  diag0.5=T
  #chosing where to save the files depending on which model you are estimating
  
  
  K_values <- list(3,4,5,6,7,8,9,10)  # Range of K values to explore
  
  print(paste0('Fitting now:' , data_description))
  
  choose_model_to_estimate = c('SST','WST','Simple')
  #-----------------------------------------------------------------------------
  # read the files in the selected folder, estimate the SST, the WST and the Simple model
  #-----------------------------------------------------------------------------
  
  for(seed in c(20,21)){
    
    
    ##############################################################################
    # Estimation: set the parameters of the estimation
    ##############################################################################
    ground_truth = NA
    n = nrow(N_ij)
    diag0.5 = F
    optimal_acceptance_rate_theta =.44
    optimal_acceptance_rate_mu = .234
    seed = seed
    N_iter <- 100000 #number of itera0tions
    burnin <- 10000 #number of discarded iterations
    thin = 5
    power_posterior_apprach = T
    if(power_posterior_apprach == F){
      saving_directory = "./Results/MCMC_output/Fixed_K/Application/"
      
      if(data_description == 'Citation_data'){
        folder_url <- "https://drive.google.com/drive/u/1/folders/1T0mrKQgDnn2QyW-NpBp2o_apFuraC9Uy"
      }else if(data_description == 'Tennis_data'){
        #where to uplY_ijoad the results
        folder_url <- "https://drive.google.com/drive/u/1/folders/1gLXPTSBCpVOZXmg9J8JxouUDHIJC1eDf"
        
      }
      
    }else if(power_posterior_apprach == T){
      saving_directory = "./Results/MCMC_output/model_choice/Powerposterior_method/"
      if(data_description == 'Citation_data'){
        folder_url <- "https://drive.google.com/drive/u/1/folders/17Zrem_-x4nD5GtLy_TmgpixVwODL57o1"
      }else if(data_description == 'Tennis_data'){
        #where to uplY_ijoad the results
        folder_url <- "https://drive.google.com/drive/u/1/folders/1VUmSgHqfLsVB1fcmqBlKAiwQHuXgMPuL"
      }
      
    }
    folder <- drive_get(as_id(folder_url))
    custom_init <- NA
    
    #number of clusters to fit
    #where to save the data
    #-----------------------------------------------------------------------------
    # SST MODEL
    #-----------------------------------------------------------------------------
    if('SST' %in% choose_model_to_estimate){
      print(paste0("Estimation of the SST model, K=",K_values))
      print(paste0("Begin cycle at:", date()))
      
      
      
      est_model = 'SST'
      
      #where to save the data
      
      
      
      #Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
      
      diag0.5<- F
      print(paste0("Estimation of the SST model, K=", K_values))
      print(paste0("Begin cycle at:", date(), "\n"))
      
      
      estimation_control <- list(z = 1, theta = 1)
      
      
      chains_SST <- adaptive_MCMC_orderstats_powerposterior(Y_ij = Y_ij, N_ij = N_ij, 
                                                            saving_directory = saving_directory,
                                                            estimation_control = estimation_control,
                                                            burnin = burnin,
                                                            ground_truth = ground_truth,
                                                            n = n, N_iter = N_iter, 
                                                            K_est = K_values,
                                                            data_description = data_description,
                                                            seed = seed, 
                                                            model = est_model, 
                                                            custom_init = custom_init,
                                                            power_posterior_apprach = power_posterior_apprach,
                                                            thin=thin,
                                                            diag0.5 = diag0.5)
      
      names(chains_SST) = paste0('chain',1:length(K_values))
      
      my_filename = paste0(saving_directory, '/Data_from',data_description, "_est_model",
                           est_model,"_Kest",paste(unlist(K_values),collapse = "_"),"seed",seed, '.rds')
      saveRDS(object = chains_SST, file = my_filename) 
      # 
      drive_put(media = my_filename, path = folder)
      #       
      #       
      #       
      #       
      
      
    }
    
    #-----------------------------------------------------------------------------
    # WST MODEL
    #-----------------------------------------------------------------------------
    
    if('WST' %in% choose_model_to_estimate){
      print(paste0("Estimation of the WST model, K=",K_values))
      print(paste0("Begin cycle at:",date()))
      #initializing each chain
      
      
      est_model = 'WST'
      
      #setting up the chain hyperparameter
      
      #Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
      
      diag0.5 <- F
      print(paste0("Estimation of the WST model, K=", K_values))
      print(paste0("Begin cycle at:", date(), "\n"))
      estimation_control <- list(z = 1,  theta = 1)
      
      chains_WST <- adaptive_MCMC_orderstats_powerposterior(Y_ij = Y_ij, N_ij = N_ij,
                                                            saving_directory = saving_directory,
                                                            estimation_control = estimation_control,
                                                            burnin = burnin,
                                                            ground_truth = ground_truth,
                                                            n = n, N_iter = N_iter, 
                                                            K_est = K_values,data_description = data_description,
                                                            seed = seed, 
                                                            model = est_model, 
                                                            custom_init = custom_init,
                                                            power_posterior_apprach = power_posterior_apprach,thin=thin,
                                                            diag0.5 = diag0.5)
      
      
      
      
      
      names(chains_WST) = paste0('chain',unlist(K_values))
      my_filename = paste0(saving_directory, '/Data_from',data_description, "_est_model",
                           est_model,"_Kest",paste(unlist(K_values),collapse = "_"),"seed",seed, '.rds')
      saveRDS(object = chains_WST, file = my_filename) 
      
      # 
      drive_put(media = my_filename, path = folder)
      # 
      
      
      
    }
    
    #-----------------------------------------------------------------------------
    # Simple model
    #-----------------------------------------------------------------------------
    
    
    if('Simple' %in% choose_model_to_estimate){
      
      print(paste0("Estimation of Simple model, K=",K_values))
      print(paste0("Begin cycle at:",date()))
      
      
      est_model = 'Simple'
      
      #setting up the chain hyperparameter
      
      
      
      
      #Boolean: power_posterior_approach = T estimates the marginal likelihood via power posteriors
      
      diag0.5 <- F
      estimation_control = list(z = 1,theta=1)
      
      chains_Simple = adaptive_MCMC_orderstats_powerposterior(Y_ij = Y_ij, N_ij = N_ij,
                                                              saving_directory = saving_directory,
                                                              estimation_control = estimation_control,
                                                              burnin = burnin,
                                                              ground_truth = ground_truth,
                                                              n = n, N_iter = N_iter, 
                                                              K_est = K_values,data_description = data_description,
                                                              seed = seed, 
                                                              model = est_model, 
                                                              custom_init = custom_init,
                                                              power_posterior_apprach = power_posterior_apprach,thin=thin,
                                                              diag0.5 = diag0.5)
      
      names(chains_Simple) = paste0('chain',1:length(K_values))
      
      my_filename = paste0(saving_directory, '/Data_from',data_description, "_est_model",
                           est_model,"_Kest",paste(unlist(K_values),collapse = "_"),"seed",seed, '.rds')
      saveRDS(object = chains_Simple, file = my_filename) 
      drive_put(media = my_filename, path = folder)
      # 
      
    }
  }
  
  
  
}