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
for(results_row in 1:nrow(results_control)){
  true_model <- ifelse(results_control$True_m[results_row], "True_ModelPOMMEst_model_","True_ModelSimpleEst_model_")
  
  filename <- list.files(pattern = paste0(true_model, est_model),path = data_wd)
  print(filename)
  
  
  
  N_iter=30000
  burnin <- 20000
  #where to save the results
  
  for(i in 1:length(filename)){
    
    #POMM model results
    est_model  <- 'POMM'
    setwd(tap)
    uploded_results<- readRDS(paste0(data_wd,"/",filename[1]))
    # Save each table with appropriate filenames and directory
    K<- nrow(uploded_results$chain1$init$P)
    K_POMM = K
    print(K)
    P<-uploded_results$chain1$ground_truth$P
    N<- nrow(uploded_results$chain1$Yij_matrix)
    N_POMM=N
    print(N)
    
    S<- if(true_model == 'True_ModelSimpleEst_model_'){
      S<-'na'
    }else{
      S<-uploded_results$chain1$ground_truth$S
    }
    z<- uploded_results$chain1$ground_truth$z
    
    P_s_title <- paste0(tap,'/P_summary_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    P_s_table_POMM <- P_summary_table(test_output = uploded_results, true_value = T, diag0.5 = TRUE, K = K, P = P, burn_in = burnin)
    mean_interval_POMM = mean(abs(P_s_table_POMM$credible_interval_95- P_s_table_POMM$credible_interval_05))
    
    
    l_op_POMM = 0 
    for(i in 1:nrow(P_s_table_POMM)){
    l_op_POMM = l_op_POMM + as.numeric(P_s_table_POMM$true_value[i] >= P_s_table_POMM$credible_interval_05[i] &&  P_s_table_POMM$true_value[i] <= P_s_table_POMM$credible_interval_95[i])
    }
    
    within_credible_POMM = l_op/nrow(P_s_table_POMM)
    mean_percentage_error_POMM = mean(abs(P_s_table_POMM$mean_est- P_s_table_POMM$true_value)/P_s_table_POMM$true_value)
    
    
    
    z_s_title <- paste0(tap,'/z_summary_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    z_s_table_POMM <- z_summary_table(test_output = uploded_results, true_value = T, diag0.5 = TRUE, K = K, z = z, burn_in = burnin)$table
    
    
    
    
    if(est_model == 'POMM'){
      S_s_title <- paste0(tap,'/S_summary_table',true_model,est_model,'_K', K,'_N', N, '.tex')
      S_s_table<- S_summary_table(test_output = uploded_results, true_value = ifelse(true_model == 'True_ModelSimpleEst_model_',F,T), diag0.5 = TRUE, S = S, K = K, burn_in = burnin)
      save_table_to_file(S_s_table, S_s_title)
    }
    
    P_d_title <- paste0(tap,'/P_diagnostic_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    P_d_table_POMM<- P_diagnostic_table(chains = uploded_results, true_value = T, diag0.5 = TRUE,K = K, P = P, burn_in = burnin,N_iter = N_iter)
    
    
    
    z_d_title <- paste0(tap,'/z_diagnostic_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    z_d_table_POMM <- z_diagnostic_table(chains = uploded_results, true_value = T, diag0.5 = TRUE, K = K, z = z, burn_in = burnin, N_iter)
    
    
    #Simple model tables
    
    est_model  <- 'Simple'
    filename <- list.files(pattern = paste0(true_model, est_model),path = data_wd)
    print(filename)
    #where to save the results
    setwd(tap)
    uploded_results<- readRDS(paste0(data_wd,"/",filename[1]))
    # Save each table with appropriate filenames and directory
    K_Simple<- nrow(uploded_results$chain1$init$P)
    K<- nrow(uploded_results$chain1$init$P)
    print(K)
    P<-uploded_results$chain1$ground_truth$P
    N_Simple<- nrow(uploded_results$chain1$Yij_matrix)
    N <- nrow(uploded_results$chain1$Yij_matrix)
    stopifnot(N_Simple == N_POMM && K_Simple == K_POMM)
    print(N)
    
    z<- uploded_results$chain1$ground_truth$z
    
    P_s_title <- paste0(tap,'/P_summary_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    P_s_table_Simple<- P_summary_table(test_output = uploded_results, true_value = T, diag0.5 = TRUE, K = K, P = P, burn_in = burnin)
    
    l_op_Simple = 0 
    for(i in 1:nrow(P_s_table_Simple)){
      l_op_Simple = l_op_Simple + as.numeric(P_s_table_Simple$true_value[i] >= P_s_table_Simple$credible_interval_05[i] &&  P_s_table_Simple$true_value[i] <= P_s_table_Simple$credible_interval_95[i])
    }
    
    within_credible_Simple = l_op_Simple/nrow(P_s_table_Simple)
    mean_percentage_error_Simple = mean(abs(P_s_table_Simple$mean_est- P_s_table_Simple$true_value)/P_s_table_Simple$true_value)
    
    mean_interval_Simple = mean(abs(P_s_table_Simple$credible_interval_95- P_s_table_Simple$credible_interval_05))
    
    z_s_title <- paste0(tap,'/z_summary_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    z_s_table_Simple<- z_summary_table(test_output = uploded_results, true_value = T, diag0.5 = TRUE, K = K, z = z, burn_in = burnin)$table
    
    
    
    P_d_title <- paste0(tap,'/P_diagnostic_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    P_d_table_Simple<- P_diagnostic_table(chains = uploded_results, true_value = T, diag0.5 = TRUE,K = K, P = P, burn_in = burnin,N_iter = N_iter)
    
    
    z_d_title <- paste0(tap,'/z_diagnostic_table',true_model,est_model,'_K', K,'_N', N, '.tex')
    z_d_table_Simple <- z_diagnostic_table(chains = uploded_results, true_value = T, diag0.5 = TRUE, K = K, z = z, burn_in = burnin, N_iter)
    
    
  }
  
  K=3
  N=100
  P_comparison <- data.frame(Model = c('POMM', 'Simple'),mean_perc_error = c(mean_percentage_error_POMM,mean_percentage_error_Simple),
             within_credible= c(within_credible_POMM,within_credible_Simple),mean_interval_length = c(mean_interval_POMM,mean_interval_Simple)
             )
  
  
  title_text <- paste0("Comparison_of_P_",'K',K,'_N',N)
  table <- gt(P_comparison, rowname_col = "Model") %>%
    tab_header(title = md(title_text)) %>%
    fmt_percent(columns = c("mean_perc_error", "within_credible")) %>%
    tab_style(
      style = cell_text(weight = "normal"),
      locations = cells_body(columns = c("mean_perc_error", "within_credible"))
    )%>%
  as_latex()%>% as.character()%>% writeLines(con = title_text)
  
  
  
  
  #fixing the label switching
  runPOMM<- label.switching(method = 'ECR' ,zpivot = obj_POMM$z_true,z = t(obj_POMM$z_container), K = K)
  # apply the permutations returned by typing:
  perm.POMM<- permute_array(array_samples = obj_POMM$p_container, perm_matrix = runPOMM$permutations$ECR)
  
  
  
  
  
  
  
  
  
  
  
  
  
  