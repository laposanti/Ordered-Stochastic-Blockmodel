


library(igraph)
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

source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/MCMC_wrapper_ORDERED.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/MCMC_wrapper_UNORDERED.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/order_statistics_model/MCMC_functions.R")




name_of_models = c('SST','WST', "Simple")
K_estimated = c(3,4,5)

#Uploading data


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

# ###############################################################################
# estimating the parameters 
# ###############################################################################

#where the data are stored
data_wd<- "/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/estimates_file/"

processed_wd <- "/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/estimates/"


#est_model  <- ifelse(results_row, "POMM","Simple")

true_model = 'Tennis_data'

for(est_model in c("SST","WST", "Simple")){
  
  filenames <- list.files(pattern = paste0('True_Model',true_model,'Est_model_', est_model),path = data_wd)
  print(filenames)
  
  for(file in 1:length(filenames)){
    
    uploaded_results<- readRDS(paste0(data_wd,"/",filenames[file]))
    print(paste0('Now estimating ', filenames[file]))
    print(paste0(length(filenames)-file,' within the same class left '))
    
    N_iter = dim(uploaded_results$chain1$est_containers$z)[[2]]
    K = dim(uploaded_results$chain1$est_containers$P)[[1]]
    
    #-------------------------------------------------------------------------------
    # P parameter estimate
    #-------------------------------------------------------------------------------
    
    
    P_s_table <- P_summary_table(chains = uploaded_results,
                                 true_value = F,
                                 diag0.5 = TRUE,
                                 K = K, P = uploaded_results$chain1$ground_truth$P,
                                 burnin = N_iter*0.25,
                                 label_switch = T)
    
    P_s_table_save <-P_s_table$table
    
    P_s_table_sum <- P_s_table_save%>%
      summarise(
        average_credible_length = mean(abs(credible_interval_95 - credible_interval_05))
      ) %>% round(3) 
    
    P_s_table_sum = P_s_table_sum %>% mutate(model = est_model)%>% mutate(n_clust = K)
    #adjusting colnames for the current number of clusters K
    if(est_model== 'SST' & file ==1){
      Pcontainer = P_s_table_sum
    }else{
      Pcontainer =  rbind(Pcontainer,P_s_table_sum)
    }
    
    
    
    #-------------------------------------------------------------------------------
    # z parameter estimate
    #-------------------------------------------------------------------------------
    
    z_tot_table<- z_summary_table(chains  = uploaded_results, true_value = F, 
                                  diag0.5 = TRUE, K = K, burnin = N_iter*0.25,
                                  label_switch = F, tap = processed_wd)
    
    z_s_table<- z_tot_table$table 
    z_s_table = z_s_table %>% mutate(model=est_model)%>% mutate(n_clust = K)
    if(est_model=='SST'&file==1){
      z_container = z_s_table}
    else{
      z_container =  rbind(z_container,z_s_table)
    }
    
    
    
    
    if(est_model !=  'Simple'){
      
      #-------------------------------------------------------------------------------
      # sigma^2 parameter estimate
      #-------------------------------------------------------------------------------
      
      sigma_squared_s_table<- sigma_squared_summary_table(chains = uploaded_results, 
                                                          true_value = F , 
                                                          diag0.5 = TRUE, K = K, burnin = N_iter*0.25)
      
      
      sigma_squared_s_table = sigma_squared_s_table %>% mutate(model=est_model)%>% mutate(n_clust = K)
      if(est_model=='SST' & file==1){
        sigma_squared_container = sigma_squared_s_table}
      else{
        sigma_squared_container =  rbind(sigma_squared_container,sigma_squared_s_table)
      }
      
      
      #-------------------------------------------------------------------------------
      # sigma^2 parameter estimate
      #-------------------------------------------------------------------------------
      
      a_s_table<- a_summary_table(chains = uploaded_results, true_value = F, 
                                  diag0.5 = TRUE, K = K, burnin = N_iter*0.25)
      
      a_s_table = a_s_table %>% mutate(model=est_model)%>% mutate(n_clust = K)
      if(est_model=='SST' & file==1){
        a_container = a_s_table
      }
      else{
        a_container =  rbind(a_container,a_s_table)
      }
      
      #-------------------------------------------------------------------------------
      # U parameter estimate
      #-------------------------------------------------------------------------------
      # 
      # 
      U_vec_s_table<- U_vec_summary_table(chains = uploaded_results, true_value = F,
                                          diag0.5 = TRUE, K = K, burnin = N_iter*0.25)
      
      U_vec_s_table = U_vec_s_table %>% mutate(model=rep(est_model,nrow(U_vec_s_table))) %>% mutate(n_clust = rep(K,nrow(U_vec_s_table)))
      if(est_model=='SST'&file==1){
        U_vec_container = U_vec_s_table
      }
      U_vec_container =  rbind(U_vec_container,U_vec_s_table)
      
      
    }
    
    
    df_traceplot = data.frame(chain = c(rep(1,(N_iter*.75)), rep(2,(N_iter*.75)), rep(3,(N_iter*.75)),rep(4,(N_iter*.75))),
                              log_likelihood = c(uploaded_results$chain1$control_containers$A[-c(1:(N_iter*0.25))],
                                                 uploaded_results$chain2$control_containers$A[-c(1:(N_iter*0.25))],
                                                 uploaded_results$chain3$control_containers$A[-c(1:(N_iter*0.25))],
                                                 uploaded_results$chain4$control_containers$A[-c(1:(N_iter*0.25))]),
                              iterations = rep(1:(N_iter*.75),4))
    df_traceplot = df_traceplot %>% mutate(chain = factor(chain, levels = 1:4))
    
    my_sexy_traceplot<- ggplot(df_traceplot, aes(x = iterations, y = log_likelihood, color = chain, group=chain))+
      geom_line(alpha = .5)+
      labs(title = "Log likelihood for the 4 chains",
           subtitle = paste0("Number of iterations: ", N_iter," || Burnin: ", N_iter*0.25), 
           x = "Iterations",
           y = "Log likelihood")+
      theme_bw()
    traceplot_name <- paste0(processed_wd,"//traceplot",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(traceplot_name,width = 500, height = 250)
    print(my_sexy_traceplot)
    dev.off()
    
    
    #-------------------------------------------------------------------------------
    # P diagnostics 
    #-------------------------------------------------------------------------------
    
    P_d_table<- P_diagnostic_table(chains = uploaded_results, true_value = F, 
                                   diag0.5 = TRUE,K = uploaded_results$chain1$init$K, 
                                   P = uploaded_results$chain1$ground_truth$P, 
                                   burnin = N_iter*0.25, N_iter = N_iter, label_switch =T)
    
    
    P_d_table_save <- P_d_table$results %>% 
      summarise(
        mean_ESS = mean(ESS),
        mean_LAG_30 = mean(LAG_30),
        mean_acceptance_rate = mean(acceptance_rate),
        median_Gelman_rubin = median(Gelman_rubin),
      ) %>% round(3)
    
    P_d_table_save = P_d_table_save %>% mutate(model= est_model)%>% mutate(n_clust = K)
    if(est_model=='SST'&file==1){
      P_d_container = P_d_table_save
    }else{
      P_d_container =  rbind(P_d_container,P_d_table_save)
    }
    
    #-------------------------------------------------------------------------------
    # z diagnostics 
    #-------------------------------------------------------------------------------
    
    z_d_table <- z_diagnostic_table(chains = uploaded_results, true_value = F, diag0.5 = TRUE, 
                                    K = K, burnin = N_iter*0.25, N_iter=N_iter*0.25,label_switch=F)
    
    z_d_table = z_d_table %>% mutate(model= est_model) %>% mutate(n_clust = K)
    
    if(est_model=='SST'&file==1){
      z_d_container = z_d_table
    }else{
      z_d_container =  rbind(z_d_container,z_d_table)
    }
    
    
    
    
    if(est_model != 'Simple'){
      
      #-------------------------------------------------------------------------------
      # sigma^2 diagnostics 
      #-------------------------------------------------------------------------------
      
      sigma_squared_d_table <- sigma_squared_diagnostic_table(chains = uploaded_results, 
                                                              true_value = F, diag0.5 = TRUE, K = K, 
                                                              burnin = N_iter*0.25, N_iter = N_iter)
      
      sigma_squared_d_table = sigma_squared_d_table %>% mutate(model= est_model) %>% mutate(n_clust = K)
      if(est_model=='SST'&file==1){
        sigma_squared_d_container = sigma_squared_d_table
      }else{
        sigma_squared_d_container =  rbind(sigma_squared_d_container,sigma_squared_d_table)
      }
      
      
      #-------------------------------------------------------------------------------
      # a diagnostics 
      #-------------------------------------------------------------------------------
      
      a_d_table <- a_diagnostic_table(chains = uploaded_results, true_value = F, diag0.5 = TRUE, 
                                      K = K, burnin = N_iter*0.25,N_iter = N_iter)
      
      a_d_table = a_d_table %>% mutate(model= est_model)%>% mutate(n_clust = K)
      if(est_model=='SST'&file==1){
        a_d_container = a_d_table
      }else{
        a_d_container =  rbind(a_d_container,a_d_table)
      }
      
      
      #-------------------------------------------------------------------------
      # U diagnostics 
      #-------------------------------------------------------------------------
      
      U_vec_d_table <- U_vec_diagnostic_table(chains = uploaded_results, true_value = F, diag0.5 = TRUE,
                                              K = K, burnin = N_iter*0.25,N_iter = N_iter)
      
      U_vec_d_table_save = U_vec_d_table$results %>% mutate(model= est_model)%>% mutate(n_clust = K)
      
      if(est_model=='SST'&file==1){
        U_vec_d_container = U_vec_d_table_save
      }else{
        U_vec_d_container =  rbind(U_vec_d_container,U_vec_d_table_save)
      }
      
    }
    
    #---------------------------------------------------------------------------
    # Saving Plots and matrices
    #---------------------------------------------------------------------------
    P_est_title <- paste0(processed_wd,'/P_est_matrix',true_model,est_model,K, '.csv')
    P_est <- round(P_s_table$P_hat,3) %>% data.frame() 
    P_est %>% write.csv(file = P_est_title)
    
    # P_true_title <- paste0(processed_wd,'/P_true_matrix',true_model,K, '.csv')
    # P_true <- round(uploaded_results$chain1$ground_truth$P,3) %>% data.frame()
    # P_true %>% write.csv(file = P_true_title)
    
    # convergence diagnostics plot -----------------------------------------------
    
    P_list <- P_d_table$plots_list
    U_list <- U_vec_d_table$plots_list
    
    #---------------------------------------------------------------------------
    #DIAGNOSTICS PLOTS FOR P 
    #---------------------------------------------------------------------------
    
    #Gelman Rubin
    plot_name <- paste0(processed_wd,"//P_gelman_rubin_plot%03d",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    # Save the plot with the constructed file name
    png(plot_name,width = 800, height = 800)
    par(mar = c(1.5, 1.5,1.5,1.5))
    gelman.plot(P_list)
    dev.off()
    
    #Crosscorrelation
    plot_name <- paste0(processed_wd,"//P_crosscorr_plot",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(plot_name,width = 800, height = 800)
    par(mfrow = c(1,1))
    crosscorr.plot(P_list)
    dev.off()
    
    #Autocorrelation
    plot_name <- paste0(processed_wd,"//P_autocorr_plot%03d",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(plot_name,width = 800, height = 800)
    autocorr.plot(P_list[[1]])
    dev.off()
    
    #---------------------------------------------------------------------------
    #DIAGNOSTICS PLOTS FOR U
    #---------------------------------------------------------------------------
    
    if(est_model != 'Simple'){
      
      #Gelman Rubin
      plot_name <- paste0(processed_wd,"//U_gelman_rubin_plot%03d",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      # Save the plot with the constructed file name
      png(plot_name,width = 800, height = 800)
      par(mar = c(1.5, 1.5,1.5,1.5))
      gelman.plot(U_list)
      dev.off()
      
      
      
      #Crosscorrelation
      plot_name <- paste0(processed_wd,"//U_crosscorr_plot",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      png(plot_name,width = 800, height = 800)
      par(mfrow = c(1,1))
      crosscorr.plot(U_list)
      dev.off()
      
      #Autocorrelation
      plot_name <- paste0(processed_wd,"//U_autocorr_plot%03d",true_model,est_model,"K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
      png(plot_name,width = 800, height = 800)
      autocorr.plot(U_list[[1]])
      dev.off()
    }
    
    #---------------------------------------------------------------------------
    # SUMMARY PLOT PLOTS FOR z
    #---------------------------------------------------------------------------
    
    # setwd(plots_dir)
    z_plot(chains = uploaded_results , true_model= true_model,
           est_model = est_model, true_value =F , 
           diag0.5 =diag0.5 , K=K, N=nrow(uploaded_results$chain1$Y_ij), z = uploaded_results$chain1$ground_truth$z ,
           burnin =  N_iter*0.25 ,label_switch = F,tap= processed_wd)
    
    z_tot_table$my_plot
    
    #---------------------------------------------------------------------------
    # COMPARISON WITH THE EXTERNAL RANKING
    #---------------------------------------------------------------------------
    
    my_names <- read.csv("/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/rawdata/players_df.csv")
    
    plot_name <- paste0(processed_wd,"//RankvsClust_Est_model",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    # Save the plot with the constructed file name
    
    g_df =  data.frame(vertex_attr(g)) %>%
      rename(player_slug= name) %>%
      left_join(players_df, by="player_slug") %>%
      mutate(degree_pl = degree(g,mode = 'out')/degree(g,mode = 'all')) %>%
      arrange()
    est_df<- data.frame(player_slug = my_names$player_slug, est_cl = z_tot_table$memb)
    combined_df<- inner_join(g_df,est_df,by = 'player_slug')
    png(plot_name,width = 800, height = 627)
    print(rank_vs_cluster(combined_df, combined_df$est_cl,est_model = est_model))
    dev.off()
    
    #---------------------------------------------------------------------------
    # CHECKING THE HOMOGENEITY OF THE CLUSTERS: HEATMAP
    #---------------------------------------------------------------------------
    
    my_names <- read.csv("/Users/lapo_santi/Desktop/Nial/MCMC_results/applications_orderstats/tennis/rawdata/players_df.csv")
    
    plot_name <- paste0(processed_wd,"//RankvsClust_Est_model",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    # Save the plot with the constructed file name
    
    g_df =  data.frame(vertex_attr(g)) %>%
      rename(player_slug= name) %>%
      left_join(players_df, by="player_slug") %>%
      mutate(degree_pl = degree(g,mode = 'out')/degree(g,mode = 'all')) %>%
      arrange()
    A = uploaded_results$chain1$Y_ij
    N = uploaded_results$chain1$N_ij
    est_df<- data.frame(player_slug = my_names$player_slug, est_cl = z_tot_table$memb)
    combined_df<- inner_join(g_df,est_df,by = 'player_slug')
    
    
    players_considered=rownames(uploaded_results$chain1$Y_ij)
    all_possible_combinations=expand.grid(players_considered,players_considered)
    
    for(i in 1:nrow(all_possible_combinations)){
      all_possible_combinations$victories[i]<- A[all_possible_combinations$Var1[i],all_possible_combinations$Var2[i]]
    }
    for(i in 1:nrow(all_possible_combinations)){
      all_possible_combinations$games_played[i]<- N[all_possible_combinations$Var1[i],all_possible_combinations$Var2[i]]
    }
    
    all_possible_combinations = all_possible_combinations %>% mutate(percentage_won = victories/games_played)
    all_possible_combinations = all_possible_combinations %>% mutate(percentage_lost = (games_played-victories)/games_played)
    
    all_possible_combinations = all_possible_combinations %>% rename(player1= Var1, player2= Var2 )
    
    all_possible_combinations = all_possible_combinations %>% 
      inner_join(combined_df, by = c("player1" = "player_slug")) %>%
      rename_at(vars(7:(7+ncol(all_possible_combinations)-2)),function(x) paste0(x,"_player1"))
    
    all_possible_combinations = all_possible_combinations %>% 
      inner_join(combined_df, by = c("player2" = "player_slug"))
    
    combined_plot<- ggplot(all_possible_combinations, aes(x=reorder(player2,est_cl), y= reorder(player1,est_cl_player1)))+
      geom_tile(aes(fill= percentage_won), color="grey9")+
      geom_ysidecol(aes(x = degree_pl_player1, color=factor(est_cl_player1))) +
      geom_xsidecol(aes(y = 1-degree_pl, color=factor(est_cl))) +
      scale_fill_gradient(low = "white", high = "red") +
      theme_bw()+ theme(legend.direction = "vertical",
                        axis.text.x = element_blank(), 
                        axis.text.y = element_blank())+
      labs(title = 'Heatmap filled with victory percentages',
           x = paste0("Players ordered by blocks"),
           y = paste0("Playersordered by blocks"),
           fill = "% victories",
           color = "Block")
    
    
    plot_name1<- paste0(processed_wd,"//Combined_plot",est_model, "_K",K,"_N",nrow(uploaded_results$chain1$Y_ij),".png")
    png(plot_name1,width = 800, height = 594)
    print(combined_plot)
    dev.off()

    # # Initialize an empty list for plots
    # plot_list <- list()
    # K= nrow(uploaded_results$chain1$est_containers$P[,,1])
    # # Loop through the groups
    # for (group1 in 1:K) {
    #   for (group2 in 1:K) {
    #     
    #     # Filter data for the current groups
    #     subsetted_K <- all_possible_combinations %>% 
    #       filter(est_cl_player1 == group1) %>%
    #       filter(est_cl == group2)
    #     
    #     # Create the ggplot object
    #     my_plot_i <- ggplot(subsetted_K, aes(x = reorder(player2, degree_pl, decreasing = F), 
    #                                          y = reorder(player1, degree_pl_player1, decreasing = T), 
    #                                          fill = percentage_won)) +
    #       geom_tile(color = 'gray') +
    #       geom_ysidecol(aes(x = percentage_won)) +
    #       scale_fill_gradient(low = "white", high = "red") +
    #       labs(title = 'Interactions matrix', 
    #            x = paste0("Players in block", group2),
    #            y = paste0("Players in block", group1),
    #            fill = "% victories") +
    #       theme_bw() +
    #       theme(legend.direction = "vertical",
    #             axis.text.x = element_blank(), 
    #             axis.text.y = element_blank(),
    #             legend.key.size = unit(0.5, "cm"),
    #             legend.position = "none")  # Remove individual legends
    #     
    #     # Convert ggplot object to grob and store in the list
    #     plot_list[[paste0("Blocks", group2, group1)]] <- ggplotGrob(my_plot_i)
    #   }
    # }
    # 
    # # Arrange and display the plots using cowplot
    # combined_plot <- plot_grid(plotlist = plot_list, ncol = K, align = 'hv')
    # 
    # # extract the legend from one of the plots
    # my_legend <- get_legend(
    #   # create some space to the left of the legend
    #   my_plot_i+ guides(color = guide_legend(nrow = K)) +
    #     theme(legend.position = "bottom")
    # )
    # # Add a common legend to the combined plot
    # combined_plot_with_legend <- plot_grid(combined_plot, my_legend, ncol = 2, rel_widths = c(5, 1)) +theme_bw()
    # 
    # # Display the combined plot with a common legend
    # print(combined_plot_with_legend)
    
    
    
    
    png(plot_name,width = 800, height = 627)
    print(rank_vs_cluster(combined_df, combined_df$est_cl,est_model = est_model))
    dev.off()
  }
}

z_container = z_container %>% arrange(n_clust) %>% 
  relocate (n_clust) %>% relocate(where(is.character)) %>% rename(K = n_clust)

z_container %>%  write.csv(file = paste0(processed_wd,"/z_container.csv"))

Pcontainer = Pcontainer %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust) 

Pcontainer %>% write.csv(file = paste0(processed_wd,"/Pcontainer.csv"))

sigma_squared_container = sigma_squared_container %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust) 

sigma_squared_container%>% write.csv(file = paste0(processed_wd,"/sigma_squared_container.csv"))

a_container = a_container %>% arrange(n_clust) %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust)

a_container %>%  write.csv(file = paste0(processed_wd,"/a_container.csv"))

U_vec_container <- U_vec_container %>% arrange(n_clust) %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust) 

U_vec_container%>% write.csv(file = paste0(processed_wd,"/U_vec_container.csv"))



z_d_container = z_d_container %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust) 

z_d_container %>% write.csv(file = paste0(processed_wd,"/z_d_container.csv"))

P_d_container = P_d_container %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust) 

P_d_container%>% write.csv(file = paste0(processed_wd,"/P_d_container.csv"))

sigma_squared_d_container = sigma_squared_d_container %>% arrange(n_clust ) %>% 
  relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust) 

sigma_squared_d_container %>% write.csv(file = paste0(processed_wd,"/sigma_squared_d_container.csv"))

a_d_container = a_d_container %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust)

a_d_container %>%  write.csv(file = paste0(processed_wd,"/a_d_container.csv"))

U_vec_d_container = U_vec_d_container %>% arrange(n_clust ) %>% relocate (n_clust) %>% 
  relocate(where(is.character)) %>% rename(K = n_clust)

U_vec_d_container %>%  write.csv(file = paste0(processed_wd,"/U_vec_d_container.csv"))













