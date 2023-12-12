


library(ggplot2)
library(abind)
library(dplyr)
library(label.switching)
library(collpcm)
library(loo)
library(gt)
library(coda)
library(mcclust)
library(LaplacesDemon)

source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/adaptive_POMM_MCMC_function.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/Inference_orderstats.R")




name_of_models = c('SST','WST', "Simple")
K_estimated = c(3,4,5)


# ###############################################################################
# estimating the parameters 
# ###############################################################################
#where the data are stored
data_wd<- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/first_simulation_study/"

#est_model  <- ifelse(results_row, "POMM","Simple")

true_model = 'POMM'

for(model_class in c("POMM","Simple")){
  
  filenames <- list.files(pattern = paste0('True_Model',true_model,'Est_model_', model_class),path = data_wd)
  print(filenames)
  
  for(file in 1:length(filenames)){
    
    uploaded_results<- readRDS(paste0(data_wd,"/",filenames[file]))
    print(paste0('Now estimating ', filenames[file]))
    print(paste0(length(filenames)-file,' within the same class left '))
    
    if(model_class == 'POMM'){
      if(uploaded_results$chain1$ground_truth$sigma_squared==0.001){
        est_model = 'SST'
      }else if(uploaded_results$chain1$ground_truth$sigma_squared==0.01){
        est_model = 'WST'
      }
    }else{
      est_model = "Simple"
    }
    N_iter = dim(uploaded_results$chain1$est_containers$z)[[2]]
    K = dim(uploaded_results$chain1$est_containers$P)[[1]]
    #-------------------------------------------------------------------------------
    # P parameter estimate
    #-------------------------------------------------------------------------------
    
    P_s_table <- P_summary_table(chains = uploaded_results,
                                 true_value = T,
                                 diag0.5 = TRUE,
                                 K = K, P = uploaded_results$chain1$ground_truth$P,
                                 burnin = N_iter*0.25,
                                 label_switch = FALSE)
    P_s_table_save <-P_s_table$table
    
    P_s_table_sum <- P_s_table_save%>%
      summarise(
        mean_mae = mean(abs(mean_est - true_value)),
        percentage_in_interval = mean(true_value >= credible_interval_05 & true_value <= credible_interval_95) * 100 ,
        average_credible_length = mean(abs(credible_interval_95 - credible_interval_05))
      ) %>% round(3) 
    
    P_s_table_sum = P_s_table_sum %>% mutate(model = est_model)%>% mutate(n_cluster = K)
    #adjusting colnames for the current number of clusters K
    if(est_model=='SST'&K==3){
      P_container = P_s_table_sum
    }else{
      P_container =  rbind(P_container,P_s_table_sum)
    }
    View(P_container)
    
    
    #-------------------------------------------------------------------------------
    # z parameter estimate
    #-------------------------------------------------------------------------------
    
    z_tot_table<- z_summary_table(chains  = uploaded_results, true_value = T, diag0.5 = TRUE, K = K, burnin = N_iter*0.25,label_switch = F)
    
    z_s_table<- z_tot_table$table 
    z_s_table = z_s_table %>% mutate(model=est_model)%>% mutate(n_clust = K)
    if(est_model=='SST'&K==3){
      z_container = z_s_table}
    else{
      z_container =  rbind(z_container,z_s_table)
    }
    
    
    
    
    if(est_model !=  'Simple'){
      
      #-------------------------------------------------------------------------------
      # sigma^2 parameter estimate
      #-------------------------------------------------------------------------------
      
      sigma_squared_s_table<- sigma_squared_summary_table(chains = uploaded_results, 
                                                          true_value = T , 
                                                          diag0.5 = TRUE, K = K, burnin = N_iter*0.25)
      
      
      sigma_squared_s_table = sigma_squared_s_table %>% mutate(model=est_model)%>% mutate(n_clust = K)
      if(est_model=='SST' & K==3){
        sigma_squared_container = sigma_squared_s_table}
      else{
        sigma_squared_container =  rbind(sigma_squared_container,sigma_squared_s_table)
      }
      
      
      #-------------------------------------------------------------------------------
      # sigma^2 parameter estimate
      #-------------------------------------------------------------------------------
      
      a_s_table<- a_summary_table(chains = uploaded_results, true_value = true_value, 
                                  diag0.5 = TRUE, K = K, burnin = N_iter*0.25)
      
      a_s_table = a_s_table %>% mutate(model=est_model)%>% mutate(n_clust = K)
      if(est_model=='SST' & K==3){
        a_container = a_s_table
      }
      else{
        a_container =  rbind(a_container,a_s_table)
      }
      #-------------------------------------------------------------------------------
      # U parameter estimate
      #-------------------------------------------------------------------------------
      
      
      U_vec_s_table<- U_vec_summary_table(chains = uploaded_results, true_value = T,
                                          diag0.5 = TRUE, K = K, burnin = N_iter*0.25)
      
      U_vec_s_table = U_vec_s_table %>% mutate(model=rep(est_model,nrow(U_vec_s_table))) %>% mutate(n_clust = rep(K,nrow(U_vec_s_table)))
      if(est_model=='SST'&K==3){
        U_vec_container = U_vec_s_table
      }
      U_vec_container =  rbind(U_vec_container,U_vec_s_table)
      
      
    }
    #-------------------------------------------------------------------------------
    # P diagnostics 
    #-------------------------------------------------------------------------------
    
    P_d_table<- P_diagnostic_table(chains = uploaded_results, true_value = T, diag0.5 = TRUE,K = nrow(uploaded_results$chain1$ground_truth$P), 
                                   P = uploaded_results$chain1$ground_truth$P, burnin = N_iter*0.25,N_iter = N_iter, label_switch =F)
    P_d_table <- P_d_table%>% 
      summarise(
        mean_ESS = mean(ESS),
        mean_LAG_30 = mean(LAG_30),
        mean_acceptance_rate = mean(acceptance_rate),
        median_Gelman_rubin = median(Gelman_rubin),
      ) %>% round(3)
    
    P_d_table = P_d_table %>% mutate(model= est_model)%>% mutate(n_clust = K)
    if(est_model=='SST'&K==3){
      P_d_container = P_d_table
    }else{
      P_d_container =  rbind(P_d_container,P_d_table)
    }
    #-------------------------------------------------------------------------------
    # z diagnostics 
    #-------------------------------------------------------------------------------
    
    z_d_table <- z_diagnostic_table(chains = uploaded_results, true_value = T, diag0.5 = TRUE, 
                                    K = K, burnin = N_iter*0.25, N_iter=N_iter*0.25,label_switch=F)
    
    z_d_table = z_d_table %>% mutate(model= est_model) %>% mutate(n_clust = K)
    if(est_model=='SST'&K==3){
      z_d_container = z_d_table
    }else{
      z_d_container =  rbind(z_d_container,z_d_table)
    }
 
    
    
    
    if(est_model != 'Simple'){
      
      #-------------------------------------------------------------------------------
      # sigma^2 diagnostics 
      #-------------------------------------------------------------------------------
      
      sigma_squared_d_table <- sigma_squared_diagnostic_table(chains = uploaded_results, 
                                                              true_value = T, diag0.5 = TRUE, K = K, 
                                                              burnin = N_iter*0.25, N_iter = N_iter)
      
      sigma_squared_d_table = sigma_squared_d_table %>% mutate(model= est_model) %>% mutate(n_clust = K)
      if(est_model=='SST'&K==3){
        sigma_squared_d_container = sigma_squared_d_table
      }else{
        sigma_squared_d_container =  rbind(sigma_squared_d_container,sigma_squared_d_table)
      }
      
      
      #-------------------------------------------------------------------------------
      # a diagnostics 
      #-------------------------------------------------------------------------------
      
      a_d_table <- a_diagnostic_table(chains = uploaded_results, true_value = T, diag0.5 = TRUE, 
                                      K = K, burnin = N_iter*0.25,N_iter = N_iter)
      
      a_d_table = a_d_table %>% mutate(model= est_model)%>% mutate(n_clust = K)
      if(est_model=='SST'&K==3){
        a_d_container = a_d_table
      }else{
        a_d_container =  rbind(a_d_container,a_d_table)
      }

      
      #-------------------------------------------------------------------------------
      # U diagnostics 
      #-------------------------------------------------------------------------------
      
      # U_vec_d_table <- U_vec_diagnostic_table(chains = uploaded_results, true_value = T, diag0.5 = TRUE, 
      #                                         K = K, burnin = N_iter*0.25,N_iter = N_iter)
      # colnames(U_vec_d_table) <- paste0(colnames(U_vec_d_table),"K",uploaded_results$chain1$ground_truth$K)
      # U_vec_d_table = U_vec_d_table %>% mutate(model= est_model)
      # U_vec_d_container =  rbind(U_vec_d_container,U_vec_d_table,by='model')
      
    }
  }
}
P_est_title <- paste0(tap,'/P_est_matrix',filename[i], '.csv')
P_est <- round(P_s_table$P_hat,3) %>% data.frame()
save_table_to_file(P_est, P_est_title,title = 'Pestmatrix',subtitle = paste0(true_model,est_model,K,N))


# setwd(plots_dir)
z_plot(chains = uploaded_results , true_model= true_model,
       est_model = est_model, true_value =T , 
       diag0.5 =diag0.5 , K=K, N=N, z = uploaded_results$chain1$ground_truth$z ,
       burnin =  burnin ,label_switch = F,tap)

est_clusters<-z_tot_table$memb

