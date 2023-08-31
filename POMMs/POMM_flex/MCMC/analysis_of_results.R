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

#setting up containers
z_summary_results = matrix(0,2,ncol = 12)
rownames(z_summary_results)<- c("POMM","Simple")

z_diagnostic_results = matrix(0,2,15)
rownames(z_diagnostic_results)<- c("POMM","Simple")

alpha_summary_results = matrix(0,1,9)
alpha_diagnostic_results= matrix(0,1,15)

S_summary_results = matrix(0,1,9)
S_diagnostic_results = matrix(0,1,15)

P_summary_results = matrix(0, 2, 9)
rownames(P_summary_results)<- c("POMM","Simple")

P_diagnostic_results= matrix(0,2,12)
rownames(P_diagnostic_results)<- c("POMM","Simple")

#Printing all possible results



#1: POMM, 0:Simple
results_control <- data.frame(True_m=c(0,0,1,1), Est_m =c(0,1,0,1))
for(results_row in 2:nrow(results_control)){
  true_model <- ifelse(results_control$True_m[results_row], "True_ModelPOMMEst_model_","True_ModelSimpleEst_model_")
  est_model  <- ifelse(results_control$Est_m[results_row], "POMM","Simple")
  filename <- list.files(pattern = paste0(true_model, est_model),path = data_wd)
  print(filename)
  
  
  file_my<- "True_ModelSimpleEst_model_Simple__N100_K3_S0.01_alpha0.5_M13_seed123.RDS"
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
    ncol(uploded_results$chain1$est_containers$z)
    if(true_model == 'True_ModelSimpleEst_model_'){
      S<-'na'
      alpha<-'na'
    }else{
      S<-uploded_results$chain1$ground_truth$S
      alpha<-uploded_results$chain1$ground_truth$alpha
    }
    z<- uploded_results$chain1$ground_truth$z
    z0<- uploded_results$chain1$init$z
    p<-uploded_results$chain1$ground_truth$P
    p0<-uploded_results$chain1$init$P
    
    P_s_table <- P_summary_table(test_output = uploded_results, 
                                 true_value = T, 
                                 diag0.5 = TRUE, 
                                 K = K, P = P, 
                                 burn_in = burnin,
                                 label_switch = F)
    P_s_table_save <-P_s_table$table
    P_s_table_save <- P_s_table_save%>% 
      summarise(
        mean_mae = mean(abs(mean_est - true_value)),
        percentage_in_interval = mean(true_value >= credible_interval_05 & true_value <= credible_interval_95) * 100 ,
        average_credible_length = mean(abs(credible_interval_95 - credible_interval_05))
      )
    if(K==3){
      P_summary_results[est_model,c(1,4,7)]<- unlist(P_s_table_save)
    }else if(K==5){
      P_summary_results[est_model,c(2,5,8)]<- unlist(P_s_table_save)
    }else{
      P_summary_results[est_model,c(3,6,9)]<- unlist(P_s_table_save)
    }
    
    
    P_est_title <- paste0(tap,'/P_est_matrix',true_model,est_model,'_K', K,'_N', N, '.tex')
    P_est <- P_s_table$P_hat %>% data.frame()
    save_table_to_file(P_est, P_est_title,title = 'Pestmatrix',subtitle = paste0(true_model,est_model,K,N))
    
    
    z_s_table<- z_summary_table(test_output = uploded_results, true_value = T, diag0.5 = TRUE, K = K, burn_in = burnin,label_switch = F)$table
    
    if(K==3){
      z_summary_results[est_model,c(1,4,7,10)]<- unlist(z_s_table)
    }else if(K==5){
      z_summary_results[est_model,c(2,5,8,11)]<- unlist(z_s_table)
    }else{
      z_summary_results[est_model,c(3,6,9,12)]<- unlist(z_s_table)
    }
    
    
    if(est_model == 'POMM'){
      
      true_value<- ifelse(true_model == 'True_ModelSimpleEst_model_',F,T)
      
      S_s_table<- S_summary_table(test_output = uploded_results, true_value = true_value , diag0.5 = TRUE, S = S, K = K, burn_in = burnin)
      
      if(true_value==T){
      if(K==3){
        S_summary_results[1,c(1,4,7)]<- unlist(S_s_table)
      }else if(K==5){
        S_summary_results[1,c(2,5,8)]<- unlist(S_s_table)
      }else{
        S_summary_results[1,c(3,6,9)]<- unlsit(S_s_table)
      }}else{
        if(K==3){
          S_summary_results[1,c(1,4)]<- unlist(S_s_table)
        }else if(K==5){
          S_summary_results[1,c(2,5)]<- unlist(S_s_table)
        }else{
          S_summary_results[1,c(3,6)]<- unlist(S_s_table)
        }
      }
      
      alpha_s_table<- alpha_summary_table(test_output = uploded_results, true_value = true_value, diag0.5 = TRUE, alpha = alpha, K = K, burn_in = burnin)
      if(true_value == T){
      if(K==3){
        alpha_summary_results[1,c(1,4,7)]<- unlist(alpha_s_table)
      }else if(K==5){
        alpha_summary_results[1,c(2,5,8)]<- unlist(alpha_s_table)
      }else{
        alpha_summary_results[1,c(3,6,9)]<- unlsit(alpha_s_table)
      }
      }}else{
        if(K==3){
          alpha_summary_results[1,c(1,4)]<- unlist(alpha_s_table)
        }else if(K==5){
          alpha_summary_results[1,c(2,5)]<- unlist(alpha_s_table)
        }else{
          alpha_summary_results[1,c(3,6)]<- unlsit(alpha_s_table)
        }
      }
        
    
    
    P_d_table<- P_diagnostic_table(chains = uploded_results, true_value = T, diag0.5 = TRUE,K = K, P = P, burn_in = burnin,N_iter = N_iter)
    P_d_table <- P_d_table%>% 
      summarise(
        mean_ESS = mean(ESS),
        mean_LAG_30 = mean(LAG_30),
        mean_Gelman_rubin = mean(Gelman_rubin),
        mean_acceptance_rate = mean(acceptance_rate)
      )
    
    if(K==3){
      P_diagnostic_results[est_model,c(1,4,7,10)]<- unlist(P_d_table)
    }else if(K==5){
      P_diagnostic_results[est_model,c(2,5,8,11)]<- unlist(P_d_table)
    }else{
      P_diagnostic_results[est_model,c(3,6,9,12)]<- unlist(P_d_table)
    }
    
    z_d_table <- z_diagnostic_table(chains = uploded_results, true_value = T, diag0.5 = TRUE, K = K, z = z, burn_in = burnin, N_iter)
    
    
    if(K==3){
      z_diagnostic_results[est_model,c(1,4,7,10,13)]<- unlist(z_d_table)
    }else if(K==5){
      z_diagnostic_results[est_model,c(2,5,8,11,14)]<- unlist(z_d_table)
    }else{
      z_diagnostic_results[est_model,c(3,6,9,12,15)]<- unlist(z_d_table)
    }
    
    if(est_model == 'POMM'){
      true_value = ifelse(true_model == 'True_ModelSimpleEst_model_',F,T)
      S_d_table <- S_diagnostic_table(chains = uploded_results, true_value = ifelse(true_model == 'True_ModelSimpleEst_model_',F,T), diag0.5 = TRUE, K = K, S = S, burn_in = burnin,N_iter = N_iter)
      if(true_value==T){
      if(K==3){
        S_diagnostic_results[1,c(1,4,7,10,13)]<- unlist(S_d_table)
      }else if(K==5){
        S_diagnostic_results[1,c(2,5,8,11,14)]<- unlist(S_d_table)
      }else{
        z_diagnostic_results[1,c(3,6,9,12,15)]<- unlist(S_d_table)
      }}else{
        if(K==3){
          S_diagnostic_results[1,c(1,4,7,10)]<- unlist(S_d_table)
        }else if(K==5){
          S_diagnostic_results[1,c(2,5,8,11)]<- unlist(S_d_table)
        }else{
          z_diagnostic_results[1,c(3,6,9,12)]<- unlist(S_d_table)
        }
      }
      
      alpha_d_table <- alpha_diagnostic_table(chains = uploded_results, true_value = true_value, diag0.5 = TRUE, K = K, alpha = alpha, burn_in = burnin,N_iter = N_iter)
       if(true_value==T){
      if(K==3){
        alpha_diagnostic_results[1,c(1,4,7,10,13)]<- unlist(alpha_d_table)
      }else if(K==5){
        alpha_diagnostic_results[1,c(2,5,8,11,14)]<- unlist(alpha_d_table)
      }else{
        alpha_diagnostic_results[1,c(3,6,9,12,15)]<- unlist(alpha_d_table)
      }
       }}else{
         if(K==3){
           alpha_diagnostic_results[1,c(1,4,7,10)]<- unlist(alpha_d_table)
         }else if(K==5){
           alpha_diagnostic_results[1,c(2,5,8,11)]<- unlist(alpha_d_table)
         }else{
           alpha_diagnostic_results[1,c(3,6,9,12)]<- unlist(alpha_d_table)
         }
    }
    setwd(plots_dir)
    z_plot(test_output =uploded_results , true_model= true_model,
           est_model = est_model, true_value =T , diag0.5 =diag0.5 , K=K, N=N, z = z ,burn_in =  burnin ,label_switch = F)
  }
  print(paste0('Hang on! Just ', nrow(results_control) - results_row, ' combinations to go!'))
}

#saving results
alpha_d_title <- paste0(tap,'/alpha_diagnostic_table', '.tex')
save_table_to_file(alpha_diagnostic_results, alpha_d_title,title = 'alphadiagnosticstable',)

S_d_title <- paste0(tap,'/S_diagnostic_table',true_model,est_model, '.tex')
save_table_to_file(S_diagnostic_results, S_d_title,title = 'Sdiagnosticstable')

z_d_title <- paste0(tap,'/z_diagnostic_table',true_model,est_model, '.tex')
save_table_to_file(z_diagnostic_results, z_d_title,title = 'zdiagnosticstable',)

P_d_title <- paste0(tap,'/P_diagnostic_table',true_model,est_model, '.tex')
save_table_to_file(P_diagnostic_results, P_d_title,title = 'Pdiagnosticstable',)

alpha_s_title <- paste0(tap,'/alpha_summary_table', '.tex')
save_table_to_file(alpha_summary_results, alpha_s_title,title = 'alphasummarytable',)

S_s_title <- paste0(tap,'/S_summary_table', '.tex')
save_table_to_file(S_summary_results, S_s_title,title = 'Ssummarytable',)

z_s_title <- paste0(tap,'/z_summary_table',true_model,est_model, '.tex')
save_table_to_file(z_summary_results, z_s_title,title = 'zsummarytable',)

P_s_title <- paste0(tap,'/P_summary_table', '.tex')
save_table_to_file(P_summary_results, P_s_title,title = 'Psummarytable')





#fixing the label switching
runPOMM<- label.switching(method = 'ECR' ,zpivot = obj_POMM$z_true,z = t(obj_POMM$z_container), K = K)
# apply the permutations returned by typing:
perm.POMM<- permute_array(array_samples = obj_POMM$p_container, perm_matrix = runPOMM$permutations$ECR)













