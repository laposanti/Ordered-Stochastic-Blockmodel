library(dplyr)
library(label.switching)
library(collpcm)
library(loo)
library(gt)
library(coda)
library(mcclust)
library(ggplot2)
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/adaptive_POMM_MCMC_function.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/Inference_functions.R")

#where the data are stored
data_wd<-  '/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Tennis application/results'

#where the data will be saved
tap <- '/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Tennis application/results/processed_results'
plots_dir<- '/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Tennis application/results/processed_results'


#Printing all possible results

#1: POMM, 0:Simple
est_model  <- 'POMM'
filename <- list.files(pattern = paste0('Tennis_application_Est_model_', est_model),path = data_wd)
print(filename)



N_iter=30000
burnin <- 20000
#where to save the results

for(i in 1:length(filename)){
  setwd(tap)
  uploded_results<- readRDS(paste0(data_wd,"/",filename))
  # Save each table with appropriate filenames and directory
  K<- nrow(uploded_results$chain1$init$P)
  print(K)
  
  N<- nrow(uploded_results$chain1$Yij_matrix)
  print(N)
  

  
  
  P_s_title <- paste0(tap,'/P_summary_table','Tennis_application_Est_model_',est_model,'_K', K,'_N', N, '.tex')
  P_s_table<- P_summary_table(test_output = uploded_results, true_value = F, diag0.5 = TRUE, K = K, P = P, burn_in = burnin)
  save_table_to_file(P_s_table, P_s_title)
  
  z_s_title <- paste0(tap,'/z_summary_table','Tennis_application_Est_model_',est_model,'_K', K,'_N', N, '.tex')
  z_s_table<- z_summary_table(test_output = uploded_results, true_value = F, diag0.5 = TRUE, K = K, z = z, burn_in = burnin)$table
  save_table_to_file(z_s_table, z_s_title)
  
  
  
  if(est_model == 'POMM'){
    S_s_title <- paste0(tap,'/S_summary_table','Tennis_application_Est_model_',est_model,'_K', K,'_N', N, '.tex')
    S_s_table<- S_summary_table(test_output = uploded_results, true_value = F, diag0.5 = TRUE, S = S, K = K, burn_in = burnin)
    save_table_to_file(S_s_table, S_s_title)
  }
  
  P_d_title <- paste0(tap,'/P_diagnostic_table','Tennis_application_Est_model_',est_model,'_K', K,'_N', N, '.tex')
  P_d_table<- P_diagnostic_table(chains = uploded_results, true_value = F, diag0.5 = TRUE,K = K, P = P, burn_in = burnin,N_iter = N_iter)
  save_table_to_file(P_d_table, P_d_title)
  
  z_d_title <- paste0(tap,'/z_diagnostic_table','Tennis_application_Est_model_',est_model,'_K', K,'_N', N, '.tex')
  z_d_table <- z_diagnostic_table(chains = uploded_results, true_value = F, diag0.5 = TRUE, K = K, z = z, burn_in = burnin,N_iter=N_iter)
  save_table_to_file(z_d_table, z_d_title)

  if(est_model == 'POMM'){
    S_d_title <- paste0(tap,'/S_diagnostic_table','Tennis_application_Est_model_',est_model,'_K', K,'_N', N, '.tex')
    S_d_table <- S_diagnostic_table(chains = uploded_results, true_value = F, diag0.5 = TRUE, K = K, S = S, burn_in = burnin,N_iter = N_iter)
    save_table_to_file(S_d_table, S_d_title)
  }
  setwd(plots_dir)
  z_plot(test_output =uploded_results , true_model= 'Tennis_application_Est_model_', 
         est_model = est_model, true_value =F , diag0.5 =diag0.5 , K=K, N=N, z = z ,burn_in =  burnin )
}

clust_est = z_summary_table(test_output = uploded_results, true_value = F, diag0.5 = TRUE, K = K, z = z, burn_in = burnin)$memb

setwd(plots_dir)

plot_name <- paste0("Est_block",est_model, "_K",K,"_N",N,".png")
# Save the plot with the constructed file name
png(plot_name,width = 800, height = 627)

# Create the ggplot plot with error bars and modifications
print(ggplot(g_df, aes(x = reorder(player_slug, median_rank), y = median_rank, color = factor(label_switch$clusters))) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = min_r, ymax = max_r), size = 1) +
  labs(x = "Player Name", y = "Ranking", title = paste0(est_model,"Estimated Block Membership and Ranking")) +
  scale_color_discrete(name = "Cluster") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    text = element_text(size = 10, family = "Arial"),
    plot.title = element_text(size = 13, face = "bold",margin = margin(r = 10)),
    plot.subtitle = element_text(size = 10, margin = margin(r = 10)),
    legend.text = element_text(size = 12),
    plot.margin = margin(20, 20, 20, 20)
  ))
# Close the device to save the plot
dev.off()





