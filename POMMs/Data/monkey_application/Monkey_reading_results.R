library(dplyr)
library(label.switching)
library(collpcm)
library(loo)
library(gt)
library(coda)
library(mcclust)
library(cowplot)
library(ggplot2)
library(igraph)

source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/adaptive_POMM_MCMC_function.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/Inference_functions.R")

#Uploading data
#-------------------------------------------------------------------------------
df_match <- read.csv('/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/monkey application/dominance.data.csv')
load("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/monkey application/df.ranks.Rdata")
df_rank <- get("df.ranks")  

#importing the data

head(df_match)


df_rank = df_rank[,-c(2)]

colnames(df_rank)<- c('name','rank')

monkey_df = monkey_df %>% filter(winner !=0 & loser != 0)
monkey_df = monkey_df %>% filter(Max.Agg.Aggressor == 'Supplant')
table(monkey_df$Max.Agg.Aggressor)

my_edges = monkey_df %>% select(winner, loser)
my_table <- as.matrix(table(my_edges))
e<-expand.grid(rownames(my_table),colnames(my_table))
for(i in 1:nrow(e)){
  e$weight[i]<- my_table[e$Var1[i], e$Var2[i]]
}
colnames(e)<- c('winner_slug','loser_slug','weight')
e = e %>% inner_join(df_rank, by = c('winner_slug'='name'))
e= e %>% rename(winner_rank = rank)
e = e %>% inner_join(df_rank, by = c('loser_slug'='name'))
e= e %>% rename(loser_rank = rank)
list_pl <- selecting_topx_andbottomy(e,bottom_x,top_x,bottom_range,upper_range)

A = list_pl$A_xvsy

#-------------------------------------------------------------------------------





#where the data are stored
data_wd<-  '/Users/lapo_santi/Desktop/Nial/MCMC_results/monkey application/'

#where the data will be saved
tap <- '/Users/lapo_santi/Desktop/Nial/MCMC_results/monkey application/processed_results/'
plots_dir<- '/Users/lapo_santi/Desktop/Nial/MCMC_results/monkey application/processed_results/'


#Printing all possible results
#1: POMM, 0:Simple


#Printing all possible results
controller=0
while(controller==0){
  simulated = F
  if(simulated==T){
    choose_model <- readline(prompt="True Model: POMM = T, Simple = F >")
  }
  true_model <- 'Monkey_data_Est_model_'
  
  #setting up containers ----------------
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
  
  corr_table = matrix(0,2,3)
  rownames(corr_table)<- c("POMM","Simple")
  #1: POMM, 0:Simple
  results_control <- data.frame(Est_m =c(0,1))
  for(results_row in c(F,T)){
    
    
    
    #---------------++++----------------
    #est_model  <- ifelse(results_row, "POMM","Simple")
    est_model<- ifelse(results_row, 'Simple', 'POMM')
    filename <- list.files(pattern = paste0(true_model, est_model),path = data_wd)
    print(filename)
    
    
    
    #where to save the results
    
    for(i in 1:length(filename)){
      setwd(tap)
      uploded_results<- readRDS(paste0(data_wd,"/",filename[i]))
      # Save each table with appropriate filenames and directory
      K<- nrow(uploded_results$chain1$init$P)
      print(K)
      N_iter=ncol(uploded_results$chain1$est_containers$z)
      burnin <- N_iter*0.5
      N<- nrow(uploded_results$chain1$Yij_matrix)
      print(N)
      ncol(uploded_results$chain1$est_containers$z)
      if(true_model == 'True_ModelPOMMEst_model_'){
        S<-uploded_results$chain1$ground_truth$S
        alpha<-uploded_results$chain1$ground_truth$alpha
      }else{
        S<-'na'
        alpha<-'na'
      }
      
      P_s_table <- P_summary_table(test_output = uploded_results, 
                                   true_value = simulated, 
                                   diag0.5 = TRUE, 
                                   K = K, P = NA, 
                                   burn_in = burnin,
                                   label_switch = F)
      P_s_table_save <-P_s_table$table
      
      
      
      P_est_title <- paste0(tap,'/P_est_matrix',true_model,est_model,'_K', K,'_N', N, '.csv')
      P_est <- round(P_s_table$P_hat,3) %>% data.frame()
      save_table_to_file(P_est, P_est_title,title = 'Pestmatrix',subtitle = paste0(true_model,est_model,K,N))
      
      z_tot_table<- z_summary_table(test_output = uploded_results, true_value = simulated, diag0.5 = TRUE, K = K, burn_in = burnin,label_switch = F)
      
      z_s_table<- z_tot_table$table
      
      if(K==2){
        z_summary_results[est_model,c(1,4,7,10)]<- round(unlist(z_s_table),2)
      }else if(K==3){
        z_summary_results[est_model,c(2,5,8,11)]<- round(unlist(z_s_table),2)
      }else if(K==4){
        z_summary_results[est_model,c(3,6,9,12)]<- round(unlist(z_s_table),2)
      }
      
      
      if(est_model == 'POMM'){
        
        true_value<- ifelse(true_model == 'True_ModelPOMMEst_model_',T,F)
        
        S_s_table<- S_summary_table(test_output = uploded_results, true_value = true_value , diag0.5 = TRUE, S = S, K = K, burn_in = burnin)
        
        if(true_value==T){
          if(K==2){
            S_summary_results[1,c(1,4,7)]<- unlist(S_s_table)
          }else if(K==3){
            S_summary_results[1,c(2,5,8)]<- unlist(S_s_table)
          }else if(K==4){
            S_summary_results[1,c(3,6,9)]<- unlist(S_s_table)
          }}else{
            if(K==2){
              S_summary_results[1,c(1,4)]<- unlist(S_s_table)
            }else if(K==3){
              S_summary_results[1,c(2,5)]<- unlist(S_s_table)
            }else if(K==4){
              S_summary_results[1,c(3,6)]<- unlist(S_s_table)
            }
          }
        
        alpha_s_table<- alpha_summary_table(test_output = uploded_results, true_value = true_value, diag0.5 = TRUE, alpha = alpha, K = K, burn_in = burnin)
        if(true_value == T){
          if(K==2){
            alpha_summary_results[1,c(1,4,7)]<- unlist(alpha_s_table)
          }else if(K==3){
            alpha_summary_results[1,c(2,5,8)]<- unlist(alpha_s_table)
          }else if(K==4){
            alpha_summary_results[1,c(3,6,9)]<- unlist(alpha_s_table)
          }
        }else{
          if(K==2){
            alpha_summary_results[1,c(1,4)]<- unlist(alpha_s_table)
          }else if(K==3){
            alpha_summary_results[1,c(2,5)]<- unlist(alpha_s_table)
          }else if(K==4){
            alpha_summary_results[1,c(3,6)]<- unlist(alpha_s_table)
          }
        }
        
      }
      
      P_d_table<- P_diagnostic_table(chains = uploded_results, true_value = simulated, diag0.5 = TRUE,K = K, P = P, burn_in = burnin,N_iter = N_iter,label_switch = T)
      P_d_table <- P_d_table%>% 
        summarise(
          mean_ESS = mean(ESS),
          mean_LAG_30 = mean(LAG_30),
          mean_acceptance_rate = mean(acceptance_rate),
          mean_Gelman_rubin = mean(Gelman_rubin),
        )
      
      if(K==2){
        P_diagnostic_results[est_model,c(1,4,7,10)]<- round(unlist(P_d_table),2)
      }else if(K==3){
        P_diagnostic_results[est_model,c(2,5,8,11)]<- round(unlist(P_d_table),2)
      }else{
        P_diagnostic_results[est_model,c(3,6,9,12)]<- round(unlist(P_d_table),2)
      }
      
      z_d_table <- z_diagnostic_table(chains = uploded_results, true_value = simulated, diag0.5 = TRUE, K = K, burn_in = burnin, N_iter)
      
      
      if(K==2){
        z_diagnostic_results[est_model,c(1,4,7,10)]<- round(unlist(z_d_table),2)
      }else if(K==3){
        z_diagnostic_results[est_model,c(2,5,8,11)]<- round(unlist(z_d_table),2)
      }else if(K==4){
        z_diagnostic_results[est_model,c(3,6,9,12)]<- round(unlist(z_d_table),2)
      }
      
      if(est_model == 'POMM'){
        true_value<- ifelse(true_model == 'True_ModelPOMMEst_model_',T,F)
        S_d_table <- S_diagnostic_table(chains = uploded_results, true_value = true_value, diag0.5 = TRUE, K = K, S = S, burn_in = burnin,N_iter = N_iter)
        if(true_value==T){
          if(K==2){
            S_diagnostic_results[1,c(1,4,7,10,13)]<- round(unlist(S_d_table),2)
          }else if(K==3){
            S_diagnostic_results[1,c(2,5,8,11,14)]<- round(unlist(S_d_table),2)
          }else if(K==4){
            S_diagnostic_results[1,c(3,6,9,12,15)]<- round(unlist(S_d_table),2)
          }}else{
            if(K==2){
              S_diagnostic_results[1,c(1,4,7,10)]<- round(unlist(S_d_table),2)
            }else if(K==3){
              S_diagnostic_results[1,c(2,5,8,11)]<- round(unlist(S_d_table),2)
            }else if(K==4){
              S_diagnostic_results[1,c(3,6,9,12)]<- round(unlist(S_d_table),2)
            }
          }
        
        alpha_d_table <- alpha_diagnostic_table(chains = uploded_results, true_value = simulated, diag0.5 = TRUE, K = K, burn_in = burnin,N_iter = N_iter)
        if(true_value==T){
          if(K==2){
            alpha_diagnostic_results[1,c(1,4,7,10,13)]<- round(unlist(alpha_d_table),2)
          }else if(K==3){
            alpha_diagnostic_results[1,c(2,5,8,11,14)]<- round(unlist(alpha_d_table),2)
          }else if(K==4){
            alpha_diagnostic_results[1,c(3,6,9,12,15)]<- round(unlist(alpha_d_table),2)
          }
        }else{
          if(K==2){
            alpha_diagnostic_results[1,c(1,4,7,10)]<- round(unlist(alpha_d_table),2)
          }else if(K==3){
            alpha_diagnostic_results[1,c(2,5,8,11)]<- round(unlist(alpha_d_table),2)
          }else if(K==4){
            alpha_diagnostic_results[1,c(3,6,9,12)]<- round(unlist(alpha_d_table),2)
          }
        }
      }
      
      setwd(plots_dir)
      z_plot(test_output =uploded_results , true_model= true_model,
             est_model = est_model, true_value =F , diag0.5 =diag0.5 , K=K, N=N, z = z ,burn_in =  burnin ,label_switch = F)
      
      est_clusters<-z_tot_table$memb
      
      plot_name <- paste0("RankvsClust_Est_model",est_model, "_K",K,"_N",N,".png")
      # Save the plot with the constructed file name
      
      
      
      est_df<- data.frame(player_slug = rownames(uploded_results$chain1$Yij_matrix), est_cl = z_tot_table$memb, rank=c(1:N))
      combined_df<- inner_join(list_pl$row_df,est_df,by = c('players_x'='player_slug'))
      
      
      my_plot<- ggplot(combined_df,aes(x=reorder(players_x,rank,decreasing=F),y= victories/games_played, fill= factor(est_cl)))+
        geom_col()+
        labs(title = paste0('POMM model results, K=','K'),
             subtitle = 'Teams are arranged according to the final ranking position',
             x = "Teams in Range [1-15]", 
             y = 'Probability of victory',fill="Estimated cluster") +
        theme_bw() +
        theme(legend.direction = "vertical") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
      png(plot_name,width = 800, height = 627)
      print(my_plot)
      dev.off()
      # if(K==3){
      #   corr_table[est_model,1]<- round(cor(combined_df$degree_pl,y = combined_df$est_cl),2)
      # }else if(K==4){
      #   corr_table[est_model,2]<- round(cor(combined_df$degree_pl,y = combined_df$est_cl),2)
      # }else{
      #   corr_table[est_model,3]<- round(cor(combined_df$degree_pl,y = combined_df$est_cl),2)
      # }
      
      
    }
    print(paste0('Hang on! Just ', nrow(results_control) - results_row, ' combinations to go!'))
  }
  
  #saving results
  alpha_d_title <- paste0(tap,'/alpha_diagnostic_table',true_model, '.csv')
  save_table_to_file(alpha_diagnostic_results, alpha_d_title,title = 'alphadiagnosticstable',)
  
  S_d_title <- paste0(tap,'/S_diagnostic_table',true_model,est_model, '.csv')
  save_table_to_file(S_diagnostic_results, S_d_title,title = 'Sdiagnosticstable')
  
  z_d_title <- paste0(tap,'/z_diagnostic_table',true_model, '.csv')
  save_table_to_file(z_diagnostic_results, z_d_title,title = 'zdiagnosticstable',)
  
  P_d_title <- paste0(tap,'/P_diagnostic_table',true_model, '.csv')
  save_table_to_file(P_diagnostic_results, P_d_title,title = 'Pdiagnosticstable',)
  
  alpha_s_title <- paste0(tap,'/alpha_summary_table',true_model, '.csv')
  save_table_to_file(alpha_summary_results, alpha_s_title,title = 'alphasummarytable',)
  
  S_s_title <- paste0(tap,'/S_summary_table',true_model, '.csv')
  save_table_to_file(S_summary_results, S_s_title,title = 'Ssummarytable',)
  
  z_s_title <- paste0(tap,'/z_summary_table',true_model, '.csv')
  save_table_to_file(z_summary_results, z_s_title,title = 'zsummarytable',)
  
  P_s_title <- paste0(tap,'/P_summary_table',true_model, '.csv')
  save_table_to_file(P_summary_results, P_s_title,title = 'Psummarytable')
  
  corr_s_title <- paste0(tap,'/Corr_summary_table',true_model, '.csv')
  save_table_to_file(corr_table, corr_s_title,title = 'Corrsummarytable')
  
  controller = controller+1
}

#where to save the results



