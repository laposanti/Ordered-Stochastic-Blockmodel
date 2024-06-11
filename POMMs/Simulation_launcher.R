

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


#setwd("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/")

source("./model_auxiliary_functions/Functions_priorSST.R")
source("./Metropolis_within_Gibbs_code.R")
source("./model_auxiliary_functions/MCMC_functions.R")





################################################################################
#                        Choose the directory
################################################################################

#chosing where to save the files



print('Simulation study for fixed K, for K=3,4,5,6')

is.simulation=T

#data.directory

data_directory = "./Data/Sim1_data/"
for(true_model in  c('Simple','WST')){

  filenames <- list.files(pattern = true_model,path =data_directory)
  print(filenames) #data to be estimated

  
  #
  choose_model_to_estimate = c('SST', 'WST', 'Simple')
  #-----------------------------------------------------------------------------
  # read the files in the selected folder, estimate the SST, the WST and the Simple model
  #-----------------------------------------------------------------------------
  
  for(file in 1:length(filenames)){
    data_to_be_estimated <- readRDS(paste0(data_directory,"/",filenames[file]))

    stopifnot(data_to_be_estimated$model == true_model)
    
    N_ij = data_to_be_estimated$N_ij
    n = nrow(N_ij)
    Y_ij = data_to_be_estimated$Y_ij
    ground_truth =data_to_be_estimated$ground_truth
    data_to_be_estimated$data_plot
    
    K= nrow(data_to_be_estimated$ground_truth$theta)
    
    print(paste0("True data--->", filenames[file]))
    
    ##############################################################################
    # Estimation: set the parameters of the estimation
    ##############################################################################
    
    n_chains = 4
    optimal_acceptance_rate_theta =.44
    optimal_acceptance_rate_mu = .234
    N_iter= 120000
    burnin = 80000
    chains_seeds = list(20,21,22,23)
    
    #-----------------------------------------------------------------------------
    # SST MODEL
    #-----------------------------------------------------------------------------
    if('SST' %in% choose_model_to_estimate){
      print(paste0("Estimation of the SST model, K=",K))
      print(paste0("Begin cycle at:", date()))
      
      #initializing each chain
      estimation_control = list(z = 1,sigma_squared=0, mu_vec=1,K=0,theta=1)

      
      K_chains = list(K,K,K,K)
      t_chains = rep(1,n_chains)
      chains_SST = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij , 
                                            estimation_control = estimation_control, 
                                            ground_truth = ground_truth, burnin=burnin,
                                            n = n, N_iter = N_iter,n_chains = n_chains, 
                                            optimal_acceptance_rate_theta = optimal_acceptance_rate_theta, 
                                            optimal_acceptance_rate_mu =optimal_acceptance_rate_mu,
                                            K = K_chains,
                                            seed = chains_seeds, model = 'SST', t= t_chains, custom_init = NA)
      
      
      my_names<- paste0("chain", 1:n_chains)
      names(chains_SST)<-my_names 


      filename_SST <- paste0("./results/simulation/",true_model,"_true//True_Model",filenames[file],"Est_model_SST","N", n,"_K", K,".RDS")
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
      K_chains = list(K,K,K,K)
      t_chains = rep(1,n_chains)
      estimation_control = list(z = 1,sigma_squared=1, mu_vec=1,K=0,theta=1)
      
      
      chains_WST = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij , 
                                            estimation_control = estimation_control, 
                                            ground_truth = ground_truth, 
                                            n = n, N_iter = N_iter,n_chains = n_chains, burnin = burnin,
                                            optimal_acceptance_rate_theta =optimal_acceptance_rate_theta, 
                                            optimal_acceptance_rate_mu =optimal_acceptance_rate_mu,
                                            K = K_chains,
                                            seed = chains_seeds, model = 'WST',t=t_chains, custom_init = NA)
      my_names <- paste0("chain", 1:n_chains)
      names(chains_WST)<-my_names 
      
      filename_WST <- paste0("./results/simulation/",true_model,"_true//True_Model",filenames[file],"Est_model_WST","_N", n,"_K", K,".RDS")
      saveRDS(chains_WST, file = filename_WST) #saving results
      beep("coin")
      
    }
    
    #-----------------------------------------------------------------------------
    # Simple model
    #-----------------------------------------------------------------------------
    
    
    if('Simple' %in% choose_model_to_estimate){
      print(paste0("Estimation of Simple model, K=",K))
      print(paste0("Begin cycle at:",date()))
 
      K_chains = list(K,K,K,K)

      t_chains = rep(1,n_chains)
      estimation_control = list(z = 1,sigma_squared=0, mu_vec=0,K=0,theta=1)
      
      chains_Simple = adaptive_MCMC_orderstats(Y_ij = Y_ij, N_ij = N_ij , 
                                               estimation_control = estimation_control, 
                                               ground_truth = ground_truth, 
                                               n = n, N_iter = N_iter,n_chains = n_chains, burnin = burnin,
                                               optimal_acceptance_rate_theta = optimal_acceptance_rate_theta , 
                                               optimal_acceptance_rate_mu =optimal_acceptance_rate_mu,
                                               K = K_chains,
                                               seed = chains_seeds, model = 'Simple',t=t_chains, custom_init = NA)
      my_names <- paste0("chain", 1:n_chains)
      names(chains_Simple)<- my_names 
      filename_Simple <- paste0("./results/simulation/",true_model,"_true//True_Model",filenames[file],"Est_model_Simple","_N", n,"_K", K,".RDS")
      saveRDS(chains_Simple, file = filename_Simple) #saving results
      beep("coin")
    }
  }
  

}


upper_tri_indices= which(upper.tri(chains_Simple$chain1$ground_truth$P, diag=T),arr.ind = T)


P_trace_df_post_switch <- do.call(rbind, lapply(1:(N_iter), function(j) {
  data.frame(iteration = j,
             P = upper.tri.extractor(chains_Simple$chain1$est_containers$theta[,,j]), 
             P_true = upper.tri.extractor(chains_Simple$chain1$ground_truth$P), 
             P_ij = paste0(upper_tri_indices[,1], upper_tri_indices[,2]))
}))


 

