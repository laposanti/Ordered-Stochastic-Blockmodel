

source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/adaptive_POMM_MCMC_function.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/Inference_functions.R")

library(EnvStats)
library(ggplot2)
library(dplyr)
library(truncnorm)
library(fossil)

library(label.switching)
library(igraph)
library(ggraph)








match_2017_url <- 'https://pkgstore.datahub.io/sports-data/atp-world-tour-tennis-data/match_scores_2017_unindexed_csv/data/df00561878fee97bf28b92cc70ae1d54/match_scores_2017_unindexed_csv.csv'
ranking_url <- 'https://pkgstore.datahub.io/sports-data/atp-world-tour-tennis-data/rankings_1973-2017_csv/data/79dd58b82401b1e872c23d4d2b6365fb/rankings_1973-2017_csv.csv'

#importing the data

df_rank <- read.csv('/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Tennis application/rankings_1973-2017_csv.csv')

match_2017_url <- 'https://pkgstore.datahub.io/sports-data/atp-world-tour-tennis-data/match_scores_2017_unindexed_csv/data/df00561878fee97bf28b92cc70ae1d54/match_scores_2017_unindexed_csv.csv'
df_match <- read.csv(match_2017_url)

head(df_match)

#computing the median rank for each player in 2017
ranks= df_rank   %>%
  filter(week_year==2017)  %>% group_by(player_slug) %>% summarise(median_rank = median(rank_number))

top100players = ranks %>% filter(median_rank <= 100) %>% arrange(median_rank)

#adding one extra column with the player id
df_r =  inner_join(ranks,df_rank%>% select(player_slug,player_id), by='player_slug')

#now, for each game I want to filter just those players in the top one-hundred
df_match = df_match %>% mutate(player_id = winner_player_id)
df_match = left_join(df_match, df_r, by = "player_id", suffix = c("_winner","_winner"))
names(df_match)[names(df_match) == "rank_number"] <- "winner_rank_number"
names(df_match)[names(df_match) == "player_slug"] <- "winner_player_slug"

df_match = df_match %>% mutate(player_id = loser_player_id)
df_match = left_join(df_match, df_r, by = "player_id", suffix = c("_loser","_loser"))
names(df_match)[names(df_match) == "rank_number"] <- "loser_rank_number"
names(df_match)[names(df_match) == "player_slug"] <- "loser_player_slug"

df_match = df_match %>% filter(loser_rank_number <= 100) %>% filter(winner_rank_number <= 100)

my_edges = df_match %>% select(winner_slug, loser_slug)

g =graph_from_edgelist(as.matrix(my_edges),directed = T)

my_name = data.frame(player_slug=vertex_attr(g)$name)

players_df = left_join(my_name, ranks %>% select(player_slug,rank_number), by="player_slug")

A = as_adjacency_matrix(g)

A= as.matrix(A)
B = matrix(0, nrow(A), ncol(A))
B[lower.tri(B)] = A[lower.tri(A)] + t(A)[lower.tri(A)]
B[upper.tri(B)] = A[upper.tri(A)] + t(A)[upper.tri(A)]

# 
# edgelist <- data.frame(df_match$winner_name, df_match$loser_name)
# reversed_edgelist <- data.frame(df_match$loser_name,df_match$winner_name)
# 
# 
# 
# 
# 
# 
# data_clean <- data_clean %>%
#   mutate(player1_number = as.integer(factor(player1, levels = unique(c(player1, player2)))),
#          player2_number = as.integer(factor(player2, levels = unique(c(player1, player2)))))
# 
# # 
# # Create a new dataframe with player names and unique IDs
# player_data <- data.frame(player = c(data_clean$player1, data_clean$player2),
#                           unique_id = c(data_clean$player1_number, data_clean$player2_number))
# 
# 
# # Remove duplicate rows (if any)
# player_data <- unique(player_data)

# Print the resulting dataframe




# 
# n_ij_matrix = matrix(0,N,N)
# M =nrow(data_clean)
# for(i in 1:M){
#   n_ij_matrix[data_clean$player1_number[i],data_clean$player2_number[i]] = data_clean$games[i]
#   n_ij_matrix[data_clean$player2_number[i],data_clean$player1_number[i]] = data_clean$games[i]
#   }
# 
# y_ij_matrix= matrix(0,N,N)
# for(i in 1:M){
#   y_ij_matrix[data_clean$player1_number[i],data_clean$player2_number[i]] = data_clean$victories[i]
#   y_ij_matrix[data_clean$player2_number[i],data_clean$player1_number[i]] = data_clean$games[i] - data_clean$victories[i]
# }

beta_max=.75


res = MCMC_POMM_Estimation(N_iter = 20000,K_true = 3,n_ij_matrix = B, y_ij_matrix =   A, beta_max = beta_max)
res_simple = MCMC_simple_model_Estimation(N_iter = 20000,K_true = 3,n_ij_matrix = B, y_ij_matrix =   A)

z_burn = res$z[,-c(1:10000)]
label_switch = label.switching(method = "DATA-BASED",z = t(z_burn),K = 3, data=rowSums(A))
point_est2 = as.vector(label_switch$clusters)

#estimates
similarity_matrix = pr_cc(z_burn)
point_est = minVI(similarity_matrix)$cl


similarity_plot(A,point_est2,point_est2)
edsimilarity_plot(similarity_matrix,point_est2,point_est2)

g <- set_vertex_attr(g, "cluster", value = point_est2)

plot(g, 
     vertex.color = point_est2,
     vertex.label = NA,
     edge.color = 'black',
     vertex.size= rowSums(A)/2,
     edge.arrow.size = .1)



g_df =  data.frame(vertex_attr(g)) %>% 
  rename( player_slug= name) %>% 
  left_join(players_df, by="player_slug")

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




