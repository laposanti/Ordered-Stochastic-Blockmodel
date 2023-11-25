


library(EnvStats)
library(ggplot2)
library(dplyr)
library(truncnorm)
library(fossil)
library(doParallel)
library(label.switching)
library(igraph)
library(ggraph)
library(foreach)
library(progressr)
library(doFuture)


source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/adaptive_POMM_MCMC_function.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/Inference_functions.R")





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




A= list_pl$A_xvsy
#A = as_adjacency_matrix(g)


Y_ij= as.matrix(A)
N_ij = matrix(0, nrow(Y_ij), ncol(Y_ij))
N_ij[lower.tri(N_ij)] = Y_ij[lower.tri(Y_ij)] + t(Y_ij)[lower.tri(Y_ij)]
N_ij[upper.tri(N_ij)] = Y_ij[upper.tri(Y_ij)] + t(Y_ij)[upper.tri(Y_ij)]

n_chains<-4
#-----------------------------------------------------------------------------
# Estimation
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
N= nrow(Y_ij)
N_iter = 40000

for(little_k in 5:5){
  K= little_k
  beta_max=0.8
  gamma_vec = rep(1/K,K)
  diag0.5=T
  alpha=.5
  S=.01
  trunc = improper_prior5(K,beta_max ,alpha,diag0.5)
  
  #------
  #POMMM
  #------
  
  
  seed=123
  #initializing each chain
  print(paste0("Estimation of POMM model, K=",K))
  init_POMM = list()
  for(chain in 1:n_chains){
    alpha0=runif(1,0.1,3)
    
    S0=runif(1,0.1,.9)
    
    trunc=improper_prior5(K,beta_max,alpha = alpha0)
    P0_POMM= simulating_overlapping_POMM_powerlaw_norm(K,alpha0,S0,trunc,beta_max,diag0.5)
    
    z0=matrix(0,N,1)
    for(item in 1:N){
      z0[item]= sample(1:K,1)
    }
    init_POMM[[chain]] =list(z = z0,alpha=alpha0,S=S0,P=P0_POMM)
  }
  
  estimation_control = list(z = 1,alpha=1,S=0,P=1)
  ground_truth= list(z = NA,alpha=NA,S=0.01,P=NA)
  hyper_params = list(K = K,beta_max =beta_max,gamma_vec = gamma_vec,diag0.5=diag0.5)
  chains_POMM = adaptive_MCMC_POMM(Yij_matrix = Y_ij,Nij_matrix = N_ij,init = init_POMM,
                                   estimation_control = estimation_control,
                                   ground_truth = ground_truth,N = N,n_chains = n_chains,
                                   N_iter = N_iter,targ_rate = .22,
                                   hyper_params =hyper_params ,seed = seed)
  my_names <- paste0("chain", 1:n_chains)
  names(chains_POMM)<-my_names 
  
  setwd('/Users/lapo_santi/Desktop/Nial/MCMC_results/kids_application/')
  filename <- paste0("kids_data_Est_model_POMM_","_N", N,"_K", K, "_seed", seed,".RDS")
  saveRDS(chains_POMM, file = filename) #saving results
  
  #------
  #Simple 
  #------
  
  
  seed=123
  print(paste0("Estimation of Simple model, K=",K))
  init_Simple = list()
  for(chain in 1:n_chains){
    P0_Simple= matrix(.5,K,K)
    P0_Simple[upper.tri(P0_Simple)]<- runif(K*(K-1)/2,0,beta_max)
    P0_Simple[lower.tri(P0_Simple)]<- 1- P0_Simple[upper.tri(P0_Simple)]
    z0=vector()
    for(i in 1:N){
      z0=append(z0, sample(1:K,1))
    }
    init_Simple[[chain]]  =list(z = z0,P=P0_Simple)
  }
  
  
  estimation_control_Simple = list(z = 1,P=1)
  ground_truth_Simple= list(z = NA,P=NA)
  hyper_params_Simple = list(K = K,beta_max =beta_max,gamma_vec = gamma_vec,diag0.5=diag0.5)
  chains_Simple = adaptive_MCMC_simple(Yij_matrix = Y_ij,Nij_matrix = N_ij,
                                       init = init_Simple,estimation_control = estimation_control_Simple,
                                       ground_truth = ground_truth_Simple,N = N,N_iter = N_iter,n_chains = n_chains,
                                       targ_rate = .22,hyper_params =hyper_params_Simple, seed = seed)
  names(chains_Simple)<-my_names 
  setwd('/Users/lapo_santi/Desktop/Nial/MCMC_results/kids_application/')
  filename_simple <- paste0("kids_data_Est_model_Simple_","_N", N,"_K", K, "_seed", seed,".RDS")
  saveRDS(chains_Simple, file = filename_simple) #saving results
  
}


label_switch = label.switching(method = "DATA-BASED",z = t(z_burn),K = 3, data=rowSums(A))
point_est2 = as.vector(label_switch$clusters)

#estimates
similarity_matrix = pr_cc(z_burn)
point_est = minVI(similarity_matrix)$cl


similarity_plot(A,point_est2,point_est2)
similarity_plot(similarity_matrix,point_est2,point_est2)

g <- set_vertex_attr(g, "cluster", value = point_est2)

plot(g, 
     vertex.color = point_est2,
     vertex.label = NA,
     edge.color = 'black',
     vertex.size= rowSums(A)/2,
     edge.arrow.size = .1)


g_df =  data.frame(vertex_attr(g)) %>% 
  rename( player_slug= name) %>% 
  left_join(players_df, by="player_slug") %>%
  mutate(degree_pl = degree(g,mode = 'out')/degree(g,mode = 'all')) %>%
  arrange()

head(g_df)




ggplot(g_df, aes(x=rank_number, y=log(rank_number), col=as.factor(cluster)))+
  geom_point()



ggplot(g_df, aes(x = player_slug, y = cluster, color = factor(cluster))) +
  geom_point(size = 3) +
  labs(x = "Player Name", y = "Cluster") +
  scale_color_discrete(name = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cluster1_data <- subset(g_df, cluster == 1)


ggplot(g_df, aes(x = player_slug, y = cluster, color = factor(cluster))) +
  geom_point(data = g_df, aes(size = rank_number)) +
  geom_text(data = cluster1_data, aes(label = player_slug), angle = 90, hjust=-.1, vjust =1.5, size = 3) +
  scale_color_discrete(name = "Cluster") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3))


ggplot(g_df, aes(x = cluster)) +
  geom_bar(aes(fill =factor(cluster)) )+
  scale_fill_manual(values = rainbow(max(g_df$cluster))) +
  scale_x_discrete()+
  xlab("Cluster") +
  ylab("Count") +
  ggtitle("Player Distribution by Cluster")

ggplot(g_df, aes(x = rank_number)) +
  geom_bar(aes(fill =factor(cluster)) )+
  scale_fill_manual(values = rainbow(max(g_df$cluster))) +
  scale_x_discrete()+
  scale_y_discrete()+
  xlab("Cluster") +
  ylab("Count") +
  ggtitle("Player Ranking by Cluster")

#############
#SimpleModel
##############

z_burn_simple = res_simple$z[,-c(1:10000)]
label_switch_simple = label.switching(method = "DATA-BASED",z = t(z_burn_simple),K = 3, data=rowSums(A))
point_est2_simple = as.vector(label_switch_simple$clusters)

#estimates
similarity_matrix_simple = pr_cc(z_burn_simple)
point_est_simple = minVI(similarity_matrix_simple)$cl

adj.rand.index(point_est2_simple,point_est2)
plotting_graphs(A, point_est2_simple, rowSums(A))

similarity_plot(A,point_est2_simple,point_est2_simple)
similarity_plot(similarity_matrix,point_est2_simple,point_est2_simple)

g_simple <- set_vertex_attr(g, "cluster", value = point_est2_simple)

plot(g_simple, 
     vertex.color = point_est2_simple,
     vertex.label = NA,
     vertex.size = rowSums(A)/2,
     edge.color = 'black',
     edge.arrow.size = 0.05)



g_df_simple =  data.frame(vertex_attr(g_simple)) %>% 
  rename( player_slug= name) %>% 
  left_join(players_df, by="player_slug")


ggplot(g_df_simple, aes(x=rank_number, y=log(rank_number), col=as.factor(cluster)))+
  geom_point()

ggplot(g_df_simple, aes(x = player_slug, y = cluster, color = factor(cluster))) +
  geom_point(size = 3) +
  labs(x = "Player Name", y = "Cluster") +
  scale_color_discrete(name = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cluster1_data <- subset(g_df_simple, cluster == 1)


ggplot(g_df_simple, aes(x = player_slug, y = cluster, color = factor(cluster))) +
  geom_point(data = g_df_simple, aes(size = rank_number)) +
  geom_text(data = cluster1_data, aes(label = player_slug), angle = 90, hjust=-.1, vjust =1.5, size = 3) +
  scale_color_discrete(name = "Cluster") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3))


ggplot(g_df_simple, aes(x = cluster)) +
  geom_bar(aes(fill =factor(cluster)) )+
  scale_fill_manual(values = rainbow(max(g_df_simple$cluster))) +
  scale_x_discrete()+
  xlab("Cluster") +
  ylab("Count") +
  ggtitle("Player Distribution by Cluster")

ggplot(g_df_simple, aes(x = rank_number)) +
  geom_bar(aes(fill =factor(cluster)) )+
  scale_fill_manual(values = rainbow(max(g_df_simple$cluster))) +
  scale_x_discrete()+
  scale_y_discrete()+
  xlab("Cluster") +
  ylab("Count") +
  ggtitle("Player Ranking by Cluster")




