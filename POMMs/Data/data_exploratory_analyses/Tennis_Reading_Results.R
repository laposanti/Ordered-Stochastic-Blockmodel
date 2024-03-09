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
match_2017_url <- 'https://pkgstore.datahub.io/sports-data/atp-world-tour-tennis-data/match_scores_2017_unindexed_csv/data/df00561878fee97bf28b92cc70ae1d54/match_scores_2017_unindexed_csv.csv'
ranking_url <- 'https://pkgstore.datahub.io/sports-data/atp-world-tour-tennis-data/rankings_1973-2017_csv/data/79dd58b82401b1e872c23d4d2b6365fb/rankings_1973-2017_csv.csv'

#importing the data

df_rank <- read.csv('/Users/lapo_santi/Desktop/Nial/raw_tennis_data/rankings_1973-2017_csv.csv')

match_2017_url <- 'https://pkgstore.datahub.io/sports-data/atp-world-tour-tennis-data/match_scores_2017_unindexed_csv/data/df00561878fee97bf28b92cc70ae1d54/match_scores_2017_unindexed_csv.csv'
df_match <- read.csv(match_2017_url)

head(df_match)

#computing the median rank for each player in 2017
ranks= df_rank   %>%
  filter(week_year==2017)  %>% group_by(player_slug) %>% summarise(median_rank = median(rank_number),max_r = max(rank_number),min_r = min(rank_number))

top100players = ranks %>% filter(median_rank <= 100) %>% arrange(median_rank)

#adding one extra column with the player id
df_r =  inner_join(ranks,df_rank%>% select(player_slug,player_id), by='player_slug')

#now, for each game I want to filter just those players in the top one-hundred
df_match = df_match %>% filter(winner_slug %in% top100players$player_slug) %>% filter(loser_slug %in% top100players$player_slug)


my_edges = df_match %>% select(winner_slug, loser_slug)
g =graph_from_edgelist(as.matrix(my_edges),directed = T)
my_name = data.frame(player_slug=vertex_attr(g)$name)

players_df = inner_join(my_name,top100players, by="player_slug")

A = as_adjacency_matrix(g)

#-------------------------------------------------------------------------------

# Assuming players_df is a data frame resulting from inner join
players_df <- inner_join(my_name, top100players, by = "player_slug")

# Save the data frame to a CSV file
write.csv(players_df, "/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/rawdata/players_df.csv", row.names = FALSE)

# Assuming A is an adjacency matrix
A <- as.matrix(as_adjacency_matrix(g))

# Save the adjacency matrix to a text file
write.table(A, "/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/rawdata/adjacency_matrix.csv", row.names = T, col.names = T)


#where the data are stored
data_wd<-  '/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Tennis application/fixing_S'

#where the data will be saved
tap <- '/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Tennis application/fixing_S/processed_results'
plots_dir<- '/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Tennis application/fixing_S/processed_results'


#Printing all possible results
#1: POMM, 0:Simple


#Printing all possible results
controller=0
while(controller==0){
  simulated = F
  if(simulated==T){
  choose_model <- readline(prompt="True Model: POMM = T, Simple = F >")
  }
  true_model <- 'Tennis_application_Est_model_'
  
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
    est_model<- 'POMM'
    filename <- list.files(pattern = paste0('Tennis_application_Est_model_', est_model),path = data_wd)
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
      
      if(K==3){
        z_summary_results[est_model,c(1,4,7,10)]<- round(unlist(z_s_table),2)
      }else if(K==4){
        z_summary_results[est_model,c(2,5,8,11)]<- round(unlist(z_s_table),2)
      }else{
        z_summary_results[est_model,c(3,6,9,12)]<- round(unlist(z_s_table),2)
      }
      
      
      if(est_model == 'POMM'){
        
        true_value<- ifelse(true_model == 'True_ModelPOMMEst_model_',T,F)
        
        S_s_table<- S_summary_table(test_output = uploded_results, true_value = true_value , diag0.5 = TRUE, S = S, K = K, burn_in = burnin)
        
        if(true_value==T){
          if(K==3){
            S_summary_results[1,c(1,4,7)]<- unlist(S_s_table)
          }else if(K==4){
            S_summary_results[1,c(2,5,8)]<- unlist(S_s_table)
          }else{
            S_summary_results[1,c(3,6,9)]<- unlist(S_s_table)
          }}else{
            if(K==3){
              S_summary_results[1,c(1,4)]<- unlist(S_s_table)
            }else if(K==4){
              S_summary_results[1,c(2,5)]<- unlist(S_s_table)
            }else{
              S_summary_results[1,c(3,6)]<- unlist(S_s_table)
            }
          }
        
        alpha_s_table<- alpha_summary_table(test_output = uploded_results, true_value = true_value, diag0.5 = TRUE, alpha = alpha, K = K, burn_in = burnin)
        if(true_value == T){
          if(K==3){
            alpha_summary_results[1,c(1,4,7)]<- unlist(alpha_s_table)
          }else if(K==4){
            alpha_summary_results[1,c(2,5,8)]<- unlist(alpha_s_table)
          }else{
            alpha_summary_results[1,c(3,6,9)]<- unlist(alpha_s_table)
          }
        }else{
          if(K==3){
            alpha_summary_results[1,c(1,4)]<- unlist(alpha_s_table)
          }else if(K==4){
            alpha_summary_results[1,c(2,5)]<- unlist(alpha_s_table)
          }else{
            alpha_summary_results[1,c(3,6)]<- unlist(alpha_s_table)
          }
        }
        
      }
      
      P_d_table<- P_diagnostic_table(chains = uploded_results, true_value = simulated, diag0.5 = TRUE,K = K, P = P, burn_in = burnin,N_iter = N_iter)
      P_d_table <- P_d_table%>% 
        summarise(
          mean_ESS = mean(ESS),
          mean_LAG_30 = mean(LAG_30),
          mean_acceptance_rate = mean(acceptance_rate),
          mean_Gelman_rubin = mean(Gelman_rubin),
        )
      
      if(K==3){
        P_diagnostic_results[est_model,c(1,4,7,10)]<- round(unlist(P_d_table),2)
      }else if(K==4){
        P_diagnostic_results[est_model,c(2,5,8,11)]<- round(unlist(P_d_table),2)
      }else{
        P_diagnostic_results[est_model,c(3,6,9,12)]<- round(unlist(P_d_table),2)
      }
      
      z_d_table <- z_diagnostic_table(chains = uploded_results, true_value = simulated, diag0.5 = TRUE, K = K, z = z, burn_in = burnin, N_iter)
      
      
      if(K==3){
        z_diagnostic_results[est_model,c(1,4,7,10)]<- round(unlist(z_d_table),2)
      }else if(K==4){
        z_diagnostic_results[est_model,c(2,5,8,11)]<- round(unlist(z_d_table),2)
      }else{
        z_diagnostic_results[est_model,c(3,6,9,12)]<- round(unlist(z_d_table),2)
      }
      
      if(est_model == 'POMM'){
        true_value<- ifelse(true_model == 'True_ModelPOMMEst_model_',T,F)
        S_d_table <- S_diagnostic_table(chains = uploded_results, true_value = true_value, diag0.5 = TRUE, K = K, S = S, burn_in = burnin,N_iter = N_iter)
        if(true_value==T){
          if(K==3){
            S_diagnostic_results[1,c(1,4,7,10,13)]<- round(unlist(S_d_table),2)
          }else if(K==4){
            S_diagnostic_results[1,c(2,5,8,11,14)]<- round(unlist(S_d_table),2)
          }else{
            S_diagnostic_results[1,c(3,6,9,12,15)]<- round(unlist(S_d_table),2)
          }}else{
            if(K==3){
              S_diagnostic_results[1,c(1,4,7,10)]<- round(unlist(S_d_table),2)
            }else if(K==4){
              S_diagnostic_results[1,c(2,5,8,11)]<- round(unlist(S_d_table),2)
            }else{
              S_diagnostic_results[1,c(3,6,9,12)]<- round(unlist(S_d_table),2)
            }
          }
        
        alpha_d_table <- alpha_diagnostic_table(chains = uploded_results, true_value = simulated, diag0.5 = TRUE, K = K, alpha = alpha, burn_in = burnin,N_iter = N_iter)
        if(true_value==T){
          if(K==3){
            alpha_diagnostic_results[1,c(1,4,7,10,13)]<- round(unlist(alpha_d_table),2)
          }else if(K==4){
            alpha_diagnostic_results[1,c(2,5,8,11,14)]<- round(unlist(alpha_d_table),2)
          }else{
            alpha_diagnostic_results[1,c(3,6,9,12,15)]<- round(unlist(alpha_d_table),2)
          }
        }else{
          if(K==3){
            alpha_diagnostic_results[1,c(1,4,7,10)]<- round(unlist(alpha_d_table),2)
          }else if(K==4){
            alpha_diagnostic_results[1,c(2,5,8,11)]<- round(unlist(alpha_d_table),2)
          }else{
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
      
      g_df =  data.frame(vertex_attr(g)) %>% 
        rename(player_slug= name) %>% 
        left_join(players_df, by="player_slug") %>%
        mutate(degree_pl = degree(g,mode = 'out')/degree(g,mode = 'all')) %>%
        arrange()
      est_df<- data.frame(player_slug = rownames(uploded_results$chain1$Yij_matrix), est_cl = z_tot_table$memb)
      combined_df<- inner_join(g_df,est_df,by = 'player_slug')
      png(plot_name,width = 800, height = 627)
      print(rank_vs_cluster(combined_df, combined_df$est_cl,est_model = est_model))
      # Close the device to save the plot
      dev.off()
      if(K==3){
        corr_table[est_model,1]<- round(cor(combined_df$degree_pl,y = combined_df$est_cl),2)
      }else if(K==4){
        corr_table[est_model,2]<- round(cor(combined_df$degree_pl,y = combined_df$est_cl),2)
      }else{
        corr_table[est_model,3]<- round(cor(combined_df$degree_pl,y = combined_df$est_cl),2)
      }
      
      
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

for(i in 1:length(filename)){
  setwd(tap)
  uploded_results<- readRDS(paste0(data_wd,"/",filename[i]))
  # Save each table with appropriate filenames and directory
  K<- nrow(uploded_results$chain1$init$P)
  print(K)
  
  N<- nrow(uploded_results$chain1$Yij_matrix)
  print(N)
  

  
  
  P_s_title <- paste0(tap,'/P_summary_table','Tennis_application_Est_model_',est_model,'_K', K,'_N', N, '.csv')
  P_s_table<- P_summary_table(test_output = uploded_results, true_value = F, diag0.5 = TRUE, K = K, P = P, burn_in = burnin)
  save_table_to_file(P_s_table, P_s_title)
  
  z_s_title <- paste0(tap,'/z_summary_table','Tennis_application_Est_model_',est_model,'_K', K,'_N', N, '.csv')
  z_s_table_tot<- z_summary_table(test_output = uploded_results, true_value = F, diag0.5 = TRUE, K = K, z = z, burn_in = burnin)
  z_s_table<- z_s_table_tot$table
  save_table_to_file(z_s_table, z_s_title)
  
  
  
  if(est_model == 'POMM'){
    S_s_title <- paste0(tap,'/S_summary_table','Tennis_application_Est_model_',est_model,'_K', K,'_N', N, '.csv')
    S_s_table<- S_summary_table(test_output = uploded_results, true_value = F, diag0.5 = TRUE, S = S, K = K, burn_in = burnin)
    save_table_to_file(S_s_table, S_s_title)
  }
  
  P_d_title <- paste0(tap,'/P_diagnostic_table','Tennis_application_Est_model_',est_model,'_K', K,'_N', N, '.csv')
  P_d_table<- P_diagnostic_table(chains = uploded_results, true_value = F, diag0.5 = TRUE,K = K, P = P, burn_in = burnin,N_iter = N_iter)
  save_table_to_file(P_d_table, P_d_title)
  
  z_d_title <- paste0(tap,'/z_diagnostic_table','Tennis_application_Est_model_',est_model,'_K', K,'_N', N, '.csv')
  z_d_table <- z_diagnostic_table(chains = uploded_results, true_value = F, diag0.5 = TRUE, K = K, z = z, burn_in = burnin,N_iter=N_iter)
  save_table_to_file(z_d_table, z_d_title)

  if(est_model == 'POMM'){
    S_d_title <- paste0(tap,'/S_diagnostic_table','Tennis_application_Est_model_',est_model,'_K', K,'_N', N, '.csv')
    S_d_table <- S_diagnostic_table(chains = uploded_results, true_value = F, diag0.5 = TRUE, K = K, S = S, burn_in = burnin,N_iter = N_iter)
    save_table_to_file(S_d_table, S_d_title)
  }
  setwd(plots_dir)
  z_plot(test_output =uploded_results , true_model= 'Tennis_application_Est_model_', 
         est_model = est_model, true_value =F , diag0.5 =diag0.5 , K=K, N=N, z = z ,burn_in =  burnin )
  
  
}

#fixing the label switching
runPOMM<- label.switching(method = 'DATA-BASED',z = t(assembling_chains(uploded_results,burnin = 20000,parameter = 'z')), K = K,data = players_df$median_rank)



# apply the permutations returned by typing:
perm.POMM<- permute_array(array_samples = obj_POMM$p_container, perm_matrix = runPOMM$permutations$ECR)

setwd(plots_dir)

plot_name <- paste0("Est_block",est_model, "_K",K,"_N",N,".png")
# Save the plot with the constructed file name
png(plot_name,width = 800, height = 627)

g_df =  data.frame(vertex_attr(g)) %>% 
  rename(player_slug= name) %>% 
  left_join(players_df, by="player_slug") %>%
  mutate(degree_pl = degree(g,mode = 'out')/degree(g,mode = 'all')) %>%
  arrange()
est_df<- data.frame(player_slug = rownames(Y_ij), est_cl = z_tot_table$memb)
combined_df<- inner_join(g_df,est_df,by = 'player_slug')
rank_vs_cluster(combined_df, combined_df$est_cl,est_model = est_model)
# Close the device to save the plot
dev.off()

