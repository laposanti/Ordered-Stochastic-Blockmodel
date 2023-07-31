library(dplyr)
library(label.switching)
library(collpcm)
library(loo)
library(gt)
library(coda)
library(mcclust)
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/adaptive_POMM_MCMC_function.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/Inference_functions.R")

#where the data are stored
data_wd<- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results_29_07/raw_results"

#where the data will be saved
tap <- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results_29_07/tables_and_plots"
plots_dir<- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/results_29_07/tables_and_plots/plots"


#Printing all possible results

#1: POMM, 0:Simple
results_control <- data.frame(True_m=c(0,0,1,1), Est_m =c(0,1,0,1))
for(results_row in 1:nrow(printing_results_control)){
  true_model <- ifelse(results_control$True_m[results_row], "True_ModelPOMMEst_model_","True_ModelSimpleEst_model_")
  est_model  <- ifelse(results_control$Est_m[results_row], "POMM","Simple")
  filename <- list.files(pattern = paste0(true_model, est_model),path = data_wd)
  print(filename)
  
  
  N_iter=30000
  burnin <- 20000
  #where to save the results
  
  for(i in 1:length(filename)){
    setwd(tap)
    uploded_results<- readRDS(paste0(data_wd,"/",filename[i]))
    # Save each table with appropriate filenames and directory
    K<- nrow(uploded_results$chain1$init$P)
    print(K)
    P<-uploded_results$chain1$ground_truth$P
    N<- nrow(uploded_results$chain1$Yij_matrix)
    print(N)
    
    S<- if(true_model == 'True_ModelSimpleEst_model_'){
      S<-'na'
    }else{
      S<-uploded_results$chain1$ground_truth$S
    }
    z<- uploded_results$chain1$ground_truth$z
    
    P_s_title <- paste0(tap,'/P_summary_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    P_s_table<- P_summary_table(test_output = uploded_results, true_value = T, diag0.5 = TRUE, K = K, P = P, burn_in = burnin)
    save_table_to_file(P_s_table, P_s_title)
    
    z_s_title <- paste0(tap,'/z_summary_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    z_s_table<- z_summary_table(test_output = uploded_results, true_value = T, diag0.5 = TRUE, K = K, z = z, burn_in = burnin)
    save_table_to_file(z_s_table, z_s_title)
    
    
    
    S<- if(est_model == 'POMM'){
    S_s_title <- paste0(tap,'/S_summary_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    S_s_table<- S_summary_table(test_output = uploded_results, true_value = ifelse(true_model == 'True_ModelSimpleEst_model_',F,T), diag0.5 = TRUE, S = S, K = K, burn_in = burnin)
    save_table_to_file(S_s_table, S_s_title)
    }
    
    P_d_title <- paste0(tap,'/P_diagnostic_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    P_d_table<- P_diagnostic_table(chains = uploded_results, true_value = T, diag0.5 = TRUE,K = K, P = P, burn_in = burnin,N_iter = N_iter)
    save_table_to_file(P_d_table, P_d_title)
    
    # z_d_title <- paste0(tap,'/z_diagnostic_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    # z_d_table <- z_diagnostic_table(chains = uploded_results, true_value = T, diag0.5 = TRUE, K = K, z = z, burn_in = burnin)
    # save_table_to_file(z_d_table, z_d_title)
    # 
    S<- if(est_model == 'POMM'){
      S_d_title <- paste0(tap,'/S_diagnostic_table',true_model,est_model,'_K', K,'_N', N, '.tex')
      S_d_table <- S_diagnostic_table(chains = uploded_results, true_value = ifelse(true_model == 'True_ModelSimpleEst_model_',F,T), diag0.5 = TRUE, K = K, S = S, burn_in = burnin,N_iter = N_iter)
      save_table_to_file(S_d_table, S_d_title)
    }
    setwd(plots_dir)
    z_plot(test_output =uploded_results , true_model= true_model, 
           est_model = est_model, true_value =T , diag0.5 =diag0.5 , K=K, N=N, z = z ,burn_in =  burnin )
  }
}


