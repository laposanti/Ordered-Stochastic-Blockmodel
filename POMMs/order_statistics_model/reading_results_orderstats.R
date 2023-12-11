

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


#where the data are stored
data_wd<- "/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/first_simulation_study/"


#where the data will be saved
tap <- "/Users/lapo_santi/Desktop/Nial/MCMC_results/simulation_study_orderstats/"
plots_dir<- "/Users/lapo_santi/Desktop/Nial/MCMC_results/simulation_study_orderstats/"




#Printing all possible results
controller=0
while(controller==0){
  
  simulated = T
  if(simulated==T){
    choose_model <- readline(prompt="True Model: POMM = T, Simple = F >")
  }
  
  true_model <- ifelse(choose_model,'POMM',"Simple")
  
  #setting up containers ----------------
  z_summary_results = matrix(0,2,ncol = 12)
  rownames(z_summary_results)<- c("POMM","Simple")
  
  z_diagnostic_results = matrix(0,2,15)
  rownames(z_diagnostic_results)<- c("POMM","Simple")
  
  a_summary_results = matrix(0,1,9)
  a_diagnostic_results= matrix(0,1,15)
  
  sigma_squared_summary_results = matrix(0,1,9)
  sigma_squared_diagnostic_results = matrix(0,1,15)
  
  U_vec_summary_results = matrix(0,5,9)
  U_vec_diagnostic_results = matrix(0,5,15)
  
  P_summary_results = matrix(0, 2, 9)
  rownames(P_summary_results)<- c("POMM","Simple")
  
  P_diagnostic_results= matrix(0,2,12)
  rownames(P_diagnostic_results)<- c("POMM","Simple")
  
  
  
  
  
  #1: POMM, 0:Simple
  results_control <- data.frame(Est_m =c(0,1))
  for(results_row in c(F,T)){
    
    
    #---------------++++----------------
    #est_model  <- ifelse(results_row, "POMM","Simple")
    est_model<- ifelse(results_row,'Simple','POMM')
    filename <- list.files(pattern = paste0('True_Model',true_model,'Est_model_', est_model),path = data_wd)
    print(filename)
    
    
    
    #where to save the results
    
    for(i in 1:length(filename)){
      
      uploded_results<- readRDS(paste0(data_wd,"/",filename[i]))
      
      if(est_model == 'POMM'){
        # Create a data frame with the x and y values
        df <- data.frame(x = 1:20000, y = t(uploded_results$chain1$est_containers$a))
        
        # Plot using ggplot
        my_plot<-ggplot(df, aes(x = x, y = y)) +
          geom_line(group=1) +
          geom_hline(yintercept = uploded_results$chain1$ground_truth$a, col = 'red') +
          labs(title = "a parameter traceplot",
               y = "a",
               x = "iterations")
        
        # Save the plot
        plot_name <- paste0(tap,"/a_parameter_traceplot",filename[i],'.png')
        # Save the plot with the constructed file name
        png(plot_name,width = 800, height = 321)
        print(my_plot) #
        dev.off()
        
        # Create a data frame with the x and y values
        df <- data.frame(x = 1:N_iter, y = t(uploded_results$chain1$est_containers$sigma_squared))
        
        # Plot using ggplot
        my_plot<-ggplot(df, aes(x = x, y = y)) +
          geom_line(group=1) +
          geom_hline(yintercept = uploded_results$chain1$ground_truth$sigma_squared, col = 'red') +
          labs(title = paste0(expression(sigma^2),"  parameter traceplot"),
               y = expression(sigma^2),
               x = "iterations")
        
        # Save the plot
        plot_name <- paste0(tap,"/sigma_squared_parameter_traceplot",filename[i],'.png')
        # Save the plot with the constructed file name
        png(plot_name,width = 800, height = 321)
        print(my_plot) #
        dev.off()
        
        df <- data.frame(
          x = rep(1:N_iter, 6),
          y = c(
            uploded_results$chain1$est_containers$P[1, 1,],
            uploded_results$chain1$est_containers$P[1, 2,],
            uploded_results$chain1$est_containers$P[1, 3,],
            uploded_results$chain1$est_containers$P[2, 2,],
            uploded_results$chain1$est_containers$P[2, 3,],
            uploded_results$chain1$est_containers$P[3, 3,]
          ),
          matrix_element = rep(c("P[1,1]", "P[1,2]", "P[1,3]", "P[2,2]", "P[2,3]", "P[3,3]"), each = N_iter)
        )
        
        # Plot using ggplot with facets
        my_plot<- ggplot(df, aes(x = x, y = y)) +
          geom_line() +
          geom_hline(aes(yintercept = c(rep(uploded_results$chain1$ground_truth$P[1,1],N_iter),
                                        rep(uploded_results$chain1$ground_truth$P[1,2],N_iter),
                                        rep(uploded_results$chain1$ground_truth$P[1,3],N_iter),
                                        rep(uploded_results$chain1$ground_truth$P[2,2],N_iter),
                                        rep(uploded_results$chain1$ground_truth$P[2,3],N_iter),
                                        rep(uploded_results$chain1$ground_truth$P[3,3],N_iter))), col = 'red') +
          labs(title = "P matrix traceplots",
               y = expression(P[p,q]),
               x = "iterations") +
          facet_wrap(~matrix_element, scales = 'free_y')
        
        plot_name <- paste0(tap,"/P_parameter_traceplot",filename[i],'.png')
        # Save the plot with the constructed file name
        png(plot_name,width = 800, height = 321)
        print(my_plot) #
        dev.off()
        
        
      }
      # Save each table with appropriate filenames and directory
      K<- nrow(uploded_results$chain1$init$P)
      print(K)
      N_iter=ncol(uploded_results$chain1$est_containers$z)
      burnin <- N_iter*0.5
      N<- nrow(uploded_results$chain1$Y_ij)
      print(N)
      ncol(uploded_results$chain1$est_containers$z)
      if(true_model == 'POMM'){
        sigma_squared<-uploded_results$chain1$ground_truth$S
        a<-uploded_results$chain1$ground_truth$a
        
      }else{
        S<-'na'
        a<-'na'
        
      }
      
      P_true_title <- paste0(tap,'/P_true_matrix',filename[i], '.csv')
      P_true<- uploded_results$chain4$ground_truth$P
      P_true <- round(P_true,3) %>% data.frame()
      save_table_to_file(P_true, P_true_title,title = 'Ptruematrix',subtitle = paste0(true_model,est_model,K,N))
      
      P_s_table <- P_summary_table(chains = uploded_results,
                                   true_value = simulated,
                                   diag0.5 = TRUE,
                                   K = K, P = P_true,
                                   burnin = burnin,
                                   label_switch = FALSE)
      P_s_table_save <-P_s_table$table
      
      P_s_table_save <- P_s_table_save%>%
        summarise(
          mean_mae = mean(abs(mean_est - true_value)),
          percentage_in_interval = mean(true_value >= credible_interval_05 & true_value <= credible_interval_95) * 100 ,
          average_credible_length = mean(abs(credible_interval_95 - credible_interval_05))
        )
      
      if(K==3){
        P_summary_results[est_model,c(1,4,7)]<- round(unlist(P_s_table_save),2)
      }else if(K==4){
        P_summary_results[est_model,c(2,5,8)]<- round(unlist(P_s_table_save),2)
      }else{
        P_summary_results[est_model,c(3,6,9)]<- round(unlist(P_s_table_save),2)
      }
      
      P_est_title <- paste0(tap,'/P_est_matrix',filename[i], '.csv')
      P_est <- round(P_s_table$P_hat,3) %>% data.frame()
      save_table_to_file(P_est, P_est_title,title = 'Pestmatrix',subtitle = paste0(true_model,est_model,K,N))
      
      
      # z_tot_table<- z_summary_table(chains  = uploded_results, true_value = simulated, diag0.5 = TRUE, K = K, burnin = burnin,label_switch = F)
      # 
      # z_s_table<- z_tot_table$table
      # 
      # if(K==3){
      #   z_summary_results[est_model,c(1,4,7,10)]<- round(unlist(z_s_table),2)
      # }else if(K==4){
      #   z_summary_results[est_model,c(2,5,8,11)]<- round(unlist(z_s_table),2)
      # }else{
      #   z_summary_results[est_model,c(3,6,9,12)]<- round(unlist(z_s_table),2)
      # }
      
      
      if(est_model == 'POMM'){
        
        true_value<- ifelse(true_model == 'True_ModelPOMMEst_model_',T,F)
        
        sigma_squared_s_table<- sigma_squared_summary_table(chains = uploded_results, 
                                                            true_value = true_value , 
                                                            diag0.5 = TRUE, K = K, burnin = burnin)
        
        if(true_value==T){
          if(K==3){
            sigma_squared_summary_results[1,c(1,4,7)]<- unlist(sigma_squared_s_table)
          }else if(K==4){
            sigma_squared_summary_results[1,c(2,5,8)]<- unlist(sigma_squared_s_table)
          }else{
            sigma_squared_summary_results[1,c(3,6,9)]<- unlist(sigma_squared_s_table)
          }}else{
            if(K==3){
              sigma_squared_summary_results[1,c(1,4)]<- unlist(sigma_squared_s_table)
            }else if(K==4){
              sigma_squared_summary_results[1,c(2,5)]<- unlist(sigma_squared_s_table)
            }else{
              sigma_squared_summary_results[1,c(3,6)]<- unlist(sigma_squared_s_table)
            }
          }
        
        a_s_table<- a_summary_table(chains = uploded_results, true_value = true_value, 
                                    diag0.5 = TRUE, K = K, burnin = burnin)
        if(true_value == T){
          if(K==3){
            a_summary_results[1,c(1,4,7)]<- unlist(a_s_table)
          }else if(K==4){
            a_summary_results[1,c(2,5,8)]<- unlist(a_s_table)
          }else{
            a_summary_results[1,c(3,6,9)]<- unlist(a_s_table)
          }
        }else{
          if(K==3){
            a_summary_results[1,c(1,4)]<- unlist(a_s_table)
          }else if(K==4){
            a_summary_results[1,c(2,5)]<- unlist(a_s_table)
          }else{
            a_summary_results[1,c(3,6)]<- unlist(a_s_table)
          }
        }
        # U_vec_s_table<- U_vec_summary_table(chains = uploded_results, true_value = true_value,
        #                             diag0.5 = TRUE, K = K, burnin = burnin)
        # if(true_value == T){
        #   if(K==3){
        #     U_vec_summary_results[c(1:3),c(1,4,7)]<- unlist(U_vec_s_table)
        #   }else if(K==4){
        #     U_vec_summary_results[c(1:4),c(2,5,8)]<- unlist(U_vec_s_table)
        #   }else{
        #     U_vec_summary_results[c(1:5),c(3,6,9)]<- unlist(U_vec_s_table)
        #   }
        # }else{
        #   if(K==3){
        #     U_vec_summary_results[c(1:3),c(1,4)]<- unlist(U_vec_s_table)
        #   }else if(K==4){
        #     U_vec_summary_results[c(1:4),c(2,5)]<- unlist(U_vec_s_table)
        #   }else{
        #     U_vec_summary_results[c(1:5),c(3,6)]<- unlist(U_vec_s_table)
        #   }
        # }
        
      }
      
      P_d_table<- P_diagnostic_table(chains = uploded_results, true_value = simulated, diag0.5 = TRUE,K = K, 
                                     P = P_true, burnin = burnin,N_iter = N_iter, label_switch =F)
      P_d_table <- P_d_table%>% 
        summarise(
          mean_ESS = mean(ESS),
          mean_LAG_30 = mean(LAG_30),
          mean_acceptance_rate = mean(acceptance_rate),
          median_Gelman_rubin = median(Gelman_rubin),
        )
      
      if(K==3){
        P_diagnostic_results[est_model,c(1,4,7,10)]<- round(unlist(P_d_table),2)
      }else if(K==4){
        P_diagnostic_results[est_model,c(2,5,8,11)]<- round(unlist(P_d_table),2)
      }else{
        P_diagnostic_results[est_model,c(3,6,9,12)]<- round(unlist(P_d_table),2)
      }
      
      z_d_table <- z_diagnostic_table(chains = uploded_results, true_value = simulated, diag0.5 = TRUE, 
                                      K = K, burnin = burnin, N_iter=N_iter,label_switch=F)
      
      
      if(K==3){
        z_diagnostic_results[est_model,c(1,4,7,10,13)]<- round(unlist(z_d_table),2)
      }else if(K==4){
        z_diagnostic_results[est_model,c(2,5,8,11,14)]<- round(unlist(z_d_table),2)
      }else{
        z_diagnostic_results[est_model,c(3,6,9,12,15)]<- round(unlist(z_d_table),2)
      }
      
      if(est_model == 'POMM'){
        true_value<- ifelse(true_model == 'POMM',T,F)
        sigma_squared_d_table <- sigma_squared_diagnostic_table(chains = uploded_results, 
                                                                true_value = true_value, diag0.5 = TRUE, K = K, 
                                                                burnin = burnin,N_iter = N_iter)
        if(true_model=='POMM'){
          if(K==3){
            sigma_squared_diagnostic_results[1,c(1,4,7,10,13)]<- round(unlist(sigma_squared_d_table),2)
          }else if(K==4){
            sigma_squared_diagnostic_results[1,c(2,5,8,11,14)]<- round(unlist(sigma_squared_d_table),2)
          }else{
            sigma_squared_diagnostic_results[1,c(3,6,9,12,15)]<- round(unlist(sigma_squared_d_table),2)
          }
        }else{
          if(K==3){
            sigma_squared_diagnostic_results[1,c(1,4,7,10)]<- round(unlist(sigma_squared_d_table),2)
          }else if(K==4){
            sigma_squared_diagnostic_results[1,c(2,5,8,11)]<- round(unlist(sigma_squared_d_table),2)
          }else{
            sigma_squared_diagnostic_results[1,c(3,6,9,12)]<- round(unlist(sigma_squared_d_table),2)
          }
        }
        
        a_d_table <- a_diagnostic_table(chains = uploded_results, true_value = true_value, diag0.5 = TRUE, K = K, burnin = burnin,N_iter = N_iter)
        if(true_model== 'POMM'){
          if(K==3){
            a_diagnostic_results[1,c(1,4,7,10,13)]<- round(unlist(a_d_table),2)
          }else if(K==4){
            a_diagnostic_results[1,c(2,5,8,11,14)]<- round(unlist(a_d_table),2)
          }else{
            a_diagnostic_results[1,c(3,6,9,12,15)]<- round(unlist(a_d_table),2)
          }
        }else{
          if(K==3){
            a_diagnostic_results[1,c(1,4,7,10)]<- round(unlist(a_d_table),2)
          }else if(K==4){
            a_diagnostic_results[1,c(2,5,8,11)]<- round(unlist(a_d_table),2)
          }else{
            a_diagnostic_results[1,c(3,6,9,12)]<- round(unlist(a_d_table),2)
          }
        }
        # U_vec_d_table <- U_vec_diagnostic_table(chains = uploded_results, true_value = true_value, diag0.5 = TRUE, K = K, burnin = burnin,N_iter = N_iter)
        # if(true_model== 'POMM'){
        #   if(K==3){
        #     U_vec_d_table[c(1:3),c(1,4,7,10,13)]<- round(unlist(U_vec_d_table),2)
        #   }else if(K==4){
        #     U_vec_d_table[c(1:4),c(2,5,8,11,14)]<- round(unlist(U_vec_d_table),2)
        #   }else{
        #     U_vec_d_table[c(1:5),c(3,6,9,12,15)]<- round(unlist(U_vec_d_table),2)
        #   }
        # }else{
        #   if(K==3){
        #     a_diagnostic_results[,c(1,4,7,10)]<- round(unlist(U_vec_d_table),2)
        #   }else if(K==4){
        #     a_diagnostic_results[,c(2,5,8,11)]<- round(unlist(U_vec_d_table),2)
        #   }else{
        #     a_diagnostic_results[,c(3,6,9,12)]<- round(unlist(U_vec_d_table),2)
        #   }
      }
    
    
    
    # setwd(plots_dir)
    z_plot(chains = uploded_results , true_model= true_model,
           est_model = est_model, true_value =T , 
           diag0.5 =diag0.5 , K=K, N=N, z = uploded_results$chain1$ground_truth$z ,
           burnin =  burnin ,label_switch = F,tap)
    
    est_clusters<-z_tot_table$memb
    
    
  }
  print(paste0('Hang on! Just ', nrow(results_control) - results_row, ' combinations to go!'))
}

#saving results
a_d_title <- paste0(tap,'/a_diagnostic_table',true_model, est_model,'.csv')
save_table_to_file(a_diagnostic_results, a_d_title,title = 'adiagnosticstable',)

U_vec_d_title <- paste0(tap,'/U_vec_diagnostic_table',true_model, est_model,'.csv')
save_table_to_file(U_vec_diagnostic_results, U_vec_d_title,title = 'U_vecdiagnosticstable',)

sigma_squared_d_title <- paste0(tap,'/sigma_squared_diagnostic_table',true_model, est_model, '.csv')
save_table_to_file(sigma_squared_diagnostic_results, sigma_squared_d_title,title = 'Sdiagnosticstable')

z_d_title <- paste0(tap,'/z_diagnostic_table',true_model, est_model,'.csv')
save_table_to_file(z_diagnostic_results, z_d_title,title = 'zdiagnosticstable',)

P_d_title <- paste0(tap,'/P_diagnostic_table',true_model, est_model, '.csv')
save_table_to_file(P_diagnostic_results, P_d_title,title = 'Pdiagnosticstable',)

a_s_title <- paste0(tap,'/a_summary_table',true_model, est_model,'.csv')
save_table_to_file(a_summary_results, a_s_title,title = 'asummarytable',)

U_vec_s_title <- paste0(tap,'/U_vec_summary_table',true_model, est_model,'.csv')
save_table_to_file(U_vec_summary_results, U_vec_s_title,title = 'U_vecsummarytable',)

sigma_squared_s_title <- paste0(tap,'/sigma_squared_summary_table',true_model, est_model, '.csv')
save_table_to_file(sigma_squared_summary_results, sigma_squared_s_title,title = 'Ssummarytable',)

z_s_title <- paste0(tap,'/z_summary_table',true_model, est_model,'.csv')
save_table_to_file(z_summary_results, z_s_title,title = 'zsummarytable',)

P_s_title <- paste0(tap,'/P_summary_table',true_model, est_model,'.csv')
save_table_to_file(P_summary_results, P_s_title,title = 'Psummarytable')


controller = controller+1
}












