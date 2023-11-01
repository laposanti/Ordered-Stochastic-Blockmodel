
# Importing dependencies


source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/adaptive_POMM_MCMC_function.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/Inference_functions.R")
library(spectrum)
library(igraph)
library(sbm)
library(heatmaply)
library(ggraph)
library(Spectrum)
library(tidygraph)
library(ggside)
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

roger_match<-df_match %>% filter (winner_slug == 'roger-federer'|loser_slug=='roger-federer')
my_edges = df_match %>% select(winner_slug, loser_slug)

e <- my_edges %>%
  data.frame() %>%
  na.omit() %>%
  group_by(winner_slug, loser_slug, .groups = 'drop') %>%
  mutate(weight = n())

e_roger = e %>% filter(winner_slug == 'roger-federer')
g =graph_from_data_frame(as.matrix(e),directed = T)

my_name = data.frame(player_slug=vertex_attr(g)$name)
players_df = inner_join(my_name,top100players, by="player_slug")

A = as.matrix(as_adjacency_matrix(g))
g_weight =graph_from_edgelist(as.matrix(my_edges),directed = T)
A_weight<- as.matrix(as_adjacency_matrix(g_weight))

Y_ij= as.matrix(A)
N_ij = matrix(0, nrow(Y_ij), ncol(Y_ij))
N_ij[lower.tri(N_ij)] = Y_ij[lower.tri(Y_ij)] + t(Y_ij)[lower.tri(Y_ij)]
N_ij[upper.tri(N_ij)] = Y_ij[upper.tri(Y_ij)] + t(Y_ij)[upper.tri(Y_ij)]

# Clustering as network




#------------
#Poisson SBM

#model fit
fit<-sbm::estimateSimpleSBM(netMat = A_weight,model='poisson',directed = T,covariates = list(A_weight))

#connecting the membership to the player name
est_df<- data.frame(player_slug<- rownames(A_weight), label<-as.numeric(fit$memberships))
colnames(est_df)<- c('player_slug','est_cl')

# Network plot

#joining the data, inserting covariates
g_tbl<-g %>% as_tbl_graph() %>% activate(edges)%>%select(-.groups)
g_tbl<-g_tbl%>% activate(nodes) %>% inner_join(est_df, by = c('name' = 'player_slug') ) %>% 
  inner_join(players_df, by = c('name' = 'player_slug'))  %>%
  mutate(degree_pl = degree(g_tbl,mode = 'out')/degree(g_tbl,mode = 'all')) %>%
  arrange()

g_tbl_df<- g_tbl %>% activate(nodes) %>% as.data.frame()
g_tbl_df<- g_tbl_df %>% mutate(loser_rank = median_rank) %>% mutate(winner_rank = median_rank)


e = e %>% inner_join(g_tbl_df %>% select(name, winner_rank), by = c('winner_slug'='name'))
e = e %>% inner_join(g_tbl_df %>% select(name, loser_rank), by = c('loser_slug'='name'))

top10_players <- g_tbl_df%>% arrange(median_rank)%>% select(name) %>%head(10)
bottom50_players<- g_tbl_df%>% arrange(desc(median_rank))%>% select(name) %>%head(50)

A_12vs50 <- matrix(0, 10, 50)
rownames(A_12vs50) <- unlist(unlist(top10_players))
colnames(A_12vs50) <- unlist(unlist(bottom50_players))

for (i in 1:nrow(e)) {
  winner <- e$winner_slug[i]
  loser <- e$loser_slug[i]
  if (winner %in% rownames(A_12vs50) && loser %in% colnames(A_12vs50)) {
    A_12vs50[winner, loser] <- e$weight[i] 
  }
}

N_12vs50 = matrix(0,10,50)
colnames(N_12vs50) <- unlist(unlist(bottom50_players))
rownames(N_12vs50) <- unlist(unlist(top10_players))
N_12vs50[lower.tri(N_12vs50)]<- A_12vs50[lower.tri(A_12vs50)] + t(A_12vs50)[lower.tri(A_12vs50)]
N_12vs50[upper.tri(N_12vs50)]<- A_12vs50[upper.tri(A_12vs50)] + t(A_12vs50)[upper.tri(A_12vs50)]

p_12vs_50<- A_12vs50/N_12vs50
colnames(p_12vs_50) <- unlist(unlist(bottom50_players))
rownames(p_12vs_50) <- unlist(unlist(top10_players))

p_10vs_10<- p_12vs_12[-which(rownames(p_12vs_12)%in%c('rafael-nadal', 'roger-federer')),-which(colnames(p_12vs_12)%in%c('rafael-nadal', 'roger-federer'))]
marginal_x_10<- colMeans(p_12vs_50,na.rm = T)
marginal_y_50<- rowMeans(p_12vs_50,na.rm = T)
boxplot(marginal_x_10,marginal_y_50)

f





data <- expand.grid(unlist(top12_players), unlist(top12_players))
colnames(data)<- c('Players_x', 'Players_y')

data$weight <- rep(0,144)
for(i in 1:144){
  for(j in 1:144){
    data$weight[i]<- p_12vs_12[data$Players_x[i],data$Players_y[i]]
  }
}

data$marginal_x <- rep(0,144)
for(i in 1:144){
  for(j in 1:144){
    data$marginal_x[i]<- mean(p_12vs_12[data$Players_x[i],],na.rm = T)
  }
}


data$marginal_y <- rep(0,144)
for(i in 1:144){
  for(j in 1:144){
    data$marginal_y[i]<- mean(p_12vs_12[,data$Players_y[i]],na.rm = T)
  }
}


data <- data %>%
  mutate(weight = replace(weight, is.nan(weight), 0))%>%
  mutate(marginal_x = replace(marginal_x, is.nan(marginal_x), 0))%>%
  mutate(marginal_y = replace(marginal_y, is.nan(marginal_y), 0))


ggplot(data, aes(x = Players_x, y = Players_y)) +
  geom_tile(aes(fill = weight))+
  scale_fill_gradient(low = "white", high = "red") +  # Adjust the colors as needed
  geom_ysidecol(aes(x =marginal_y, fill= marginal_y))  +
  labs(title = "Heatmap with Densities", x = "Players", y = "Players") +
  theme_bw()+
  theme(legend.direction = "horizontal") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggplot(data, aes(x = Players_x, y = Players_y)) +
  geom_tile(aes(fill = weight))+
  scale_fill_gradient(low = "white", high = "red") +  # Adjust the colors as needed
  geom_ysidetile(aes(x = "Marginal y", yfill =marginal_y))  +
  geom_xsidedensity(aes(y=data %>% group_by(Players_x)%>% summarise(marginal_x = mean(marginal_x)) %>% select(marginal_x)), orientation = 'y')  +
  labs(title = "Heatmap with Densities", x = "Players", y = "Players") +
  theme_bw()+
  theme(legend.direction = "horizontal") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


marginal_x <- colMeans(p_12vs_12,na.rm = T) 
marginal_y <- rowMeans(p_12vs_12,na.rm = T) 

# Create a dataframe for plotting
plot(density(marginal_data_x$Probability))
boxplot(marginal_x,marginal_y)


marginal_data_y <- data.frame(Player = rownames(p_12vs_12), Probability = marginal_y)
plot(density(marginal_data_y$Probability))


# Plot the marginal probabilities along the x-axis
ggplot(data, aes(x = Players_x)) +
  geom_point(aes(fill = weight)) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Heatmap with Densities", x = "Players", y = "Players") +
  theme_bw() +
  theme(legend.direction = "horizontal", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

g_tbl <- g_tbl  %>%mutate(est_cl = as.factor(est_cl))


color_palette <- colorRampPalette(c("grey", "black"))(n = 5)
set.seed(123)
ggraph(g_tbl, layout = 'fr') +
  geom_edge_link(aes(alpha = weight)) + # Set custom color palette
  geom_node_point(aes(color = est_cl), size = scale(rowSums(A))+2) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank())
g_tbl_df<- g_tbl %>% activate(nodes) %>% as.data.frame()
rank_vs_cluster(g_tbl_df, g_tbl_df$est_cl,est_model = est_model)

#table data---

g_tbl_df%>%group_by(est_cl)%>% summarise(mean(degree_pl))
cor(g_tbl_df$median_rank,as.numeric(g_tbl_df$est_cl))

K_vector = vector()
corr_vector<- vector()
#computing significance of clustering
N=50
for(i in 1:N){
  A_rewired_constant_degree = matrix(0, 95,95)
  for(i in 1:95){
    
    my_max<- max(sample(A_weight[i,-i]))
    new_counts<-vector()
    while(sum(new_counts)<sum(A_weight[i,-i])){
      new_counts<- append(new_counts, sample(1:my_max),1)
      if(sum(new_counts)==sum(A_weight[i,-i])){
        break
      }else if(sum(new_counts)>sum(A_weight[i,-i])){
        new_counts<- new_counts[-length(new_counts)]
        new_counts<- append(new_counts, abs(sum(A_weight[i,-i]) - sum(new_counts)))
        break}
    }
    print(new_counts)
    random_indices<- sample(setdiff(1:95,i),length(new_counts), replace = F)
    A_rewired_constant_degree[i,random_indices] <- new_counts
  }
  
  fit_random<-sbm::estimateSimpleSBM(netMat = A_rewired_constant_degree,model='poisson')
  K<- length(unique(fit_random$memberships))
  print(K)
  my_cor<- cor(fit_random$memberships, g_tbl_df$median_rank)
  K_vector<- append(K_vector,K)
  corr_vector<- append(corr_vector,my_cor)
}
mean(K_vector)

K_vector_density = vector()
corr_vector_density<- vector()

N=50
diag_indicator<- matrix(0, 95,2)
for(j in 1:95){
  diag_indicator[j,] <- c(j,j)
}
for(i in 1:N){
  A_rewired_constant_density = matrix(0, 95,95)


  A_rewired_constant_degree[- c(diag_indicator)] <- sample(A_weight[- c(diag_indicator)], length(A_weight[- c(diag_indicator)]), replace = F)
  
  
  fit_random<-sbm::estimateSimpleSBM(netMat = A_rewired_constant_degree,model='poisson')
  K<- length(unique(fit_random$memberships))
  print(K)
  my_cor<- cor(fit_random$memberships, g_tbl_df$median_rank)
  K_vector_density<- append(K_vector_density,K)
  corr_vector_density<- append(corr_vector_density,my_cor)
}


library(factoextra)

#-------------------
#Spectral clustering

X_sc <- spectral_clustering(A)
heatmap(X_sc)

elbow_plot <- fviz_nbclust(X_sc, kmeans,method = 'wss')
silhouette_plot<- fviz_nbclust(X_sc, kmeans,method = 'sil')
n_clust_elbow<- which.max(-1*diff(elbow_plot$data$y))+1
n_clust_silhouette<- which.max(silhouette_plot$data$y)
# run kmeans on the 2 eigenvectors
X_sc_kmeans <- kmeans(X_sc, n_clust_silhouette)

#connecting the membership to the player name
est_df<- data.frame(player_slug<- rownames(A), label<-as.numeric(X_sc_kmeans$cluster))
colnames(est_df)<- c('player_slug','est_cl')

# Network plot

#joining the data, inserting covariates
g_tbl<-g %>% as_tbl_graph() %>% activate(edges)%>%select(-.groups)
g_tbl<-g_tbl%>% activate(nodes) %>% inner_join(est_df, by = c('name' = 'player_slug') ) %>% 
  inner_join(players_df, by = c('name' = 'player_slug'))  %>%
  mutate(degree_pl = degree(g_tbl,mode = 'out')/degree(g_tbl,mode = 'all')) %>%
  arrange()


g_tbl <- g_tbl %>%mutate(est_cl = as.factor(est_cl))

color_palette <- colorRampPalette(c("grey", "black"))(n = 5)
set.seed(123)
ggraph(g_tbl, layout = 'fr') +
  geom_edge_link(aes(alpha = weight)) + # Set custom color palette
  geom_node_point(aes(color = est_cl), size = scale(rowSums(A))+2) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank())
g_tbl_df<- g_tbl %>% activate(nodes) %>% as.data.frame()
rank_vs_cluster(g_tbl_df, g_tbl_df$est_cl,est_model = est_model)

#table data---

g_tbl_df%>%group_by(est_cl)%>% summarise(mean(degree_pl))
cor(g_tbl_df$median_rank,as.numeric(g_tbl_df$est_cl))

N=50
for(j in 1:N){
  A_rewired_constant_degree = matrix(0, 95,95)
  for(i in 1:95){
    my_max<- max(sample(A_weight[i,-i]))
    new_counts<-vector()
    while(sum(new_counts)<sum(A_weight[i,-i])){
      new_counts<- append(new_counts, sample(1:my_max),1)
      if(sum(new_counts)==sum(A_weight[i,-i])){
        break
      }else if(sum(new_counts)>sum(A_weight[i,-i])){
        new_counts<- new_counts[-length(new_counts)]
        new_counts<- append(new_counts, abs(sum(A_weight[i,-i]) - sum(new_counts)))
        break}
    }
    print(new_counts)
    random_indices<- sample(setdiff(1:95,i),length(new_counts), replace = F)
    A_rewired_constant_degree[i,random_indices] <- new_counts
  }
  
  X_sc_degree <- spectral_clustering(A_rewired_constant_degree)
  
  elbow_plot <- fviz_nbclust(X_sc_degree, kmeans,method = 'wss')
  silhouette_plot<- fviz_nbclust(X_sc_degree, kmeans,method = 'sil')
  n_clust_elbow<- which.max(-1*diff(elbow_plot$data$y))+1
  n_clust_silhouette<- which.max(silhouette_plot$data$y)
  # run kmeans on the 2 eigenvectors
  X_sc_kmeans <- kmeans(X_sc_degree, n_clust_silhouette)
  K<- length(unique(X_sc_kmeans$cluster))
  print(K)
  my_cor<- cor(X_sc_kmeans$cluster, g_tbl_df$median_rank)
  K_vector<- append(K_vector,K)
  corr_vector<- append(corr_vector,my_cor)
}


K_vector_density = vector()
corr_vector_density<- vector()

N=50
diag_indicator<- matrix(0, 95,2)
for(j in 1:95){
  diag_indicator[j,] <- c(j,j)
}
for(i in 1:N){
  A_rewired_constant_density = matrix(0, 95,95)
  
  
  A_rewired_constant_degree[- c(diag_indicator)] <- sample(A_weight[- c(diag_indicator)], length(A_weight[- c(diag_indicator)]), replace = F)
  
  
  X_sc_degree <- spectral_clustering(A_rewired_constant_degree)
  
  elbow_plot <- fviz_nbclust(X_sc_degree, kmeans,method = 'wss')
  silhouette_plot<- fviz_nbclust(X_sc_degree, kmeans,method = 'sil')
  n_clust_elbow<- which.max(-1*diff(elbow_plot$data$y))+1
  n_clust_silhouette<- which.max(silhouette_plot$data$y)
  # run kmeans on the 2 eigenvectors
  X_sc_kmeans <- kmeans(X_sc_degree, n_clust_silhouette)
  K<- length(unique(X_sc_kmeans$cluster))
  print(K)
  my_cor<- cor(X_sc_kmeans$cluster, g_tbl_df$median_rank)
  K_vector_density<- append(K_vector_density,K)
  corr_vector_density<- append(corr_vector_density,my_cor)
}
#===============================================================================
#Distance method
#===============================================================================
#----------------------------------------------
#Spectral clustering applied to distance matrix

d_matrix_euc<- as.matrix(dist(A_weight,upper = T,method = 'euc'))
d_matrix_man<- as.matrix(dist(A_weight,upper = T,method = 'man'))
d_matrix_jac<- as.matrix(dist(A_weight,upper = T,method = 'bin'))

dist_X_sc_euc <- spectral_clustering(d_matrix_euc)
dist_X_sc_man <- spectral_clustering(d_matrix_man)
dist_X_sc_jac <- spectral_clustering(d_matrix_jac)



silhouette_plot_euc<-fviz_nbclust(dist_X_sc_euc, kmeans,method = 'sil')
silhouette_plot_man<-fviz_nbclust(dist_X_sc_man, kmeans,method = 'sil')
silhouette_plot_jac<-fviz_nbclust(dist_X_sc_jac, kmeans,method = 'sil')

n_clust_silhouette_euc<- which.max(silhouette_plot_euc$data$y)
n_clust_silhouette_man<- which.max(silhouette_plot_man$data$y)
n_clust_silhouette_jac<- which.max(silhouette_plot_jac$data$y)

# Spectral Clustering
X_sc_kmeans_dist_euc <- kmeans(dist_X_sc_euc, n_clust_silhouette_euc)
X_sc_kmeans_dist_man<- kmeans(dist_X_sc_man, n_clust_silhouette_man)
X_sc_kmeans_dist_jac<- kmeans(dist_X_sc_jac, n_clust_silhouette_jac)

#K-means Clustering 
kmeans_dist_euc <- kmeans(d_matrix_euc, n_clust_silhouette_euc)
kmeans_dist_man<- kmeans(d_matrix_man, n_clust_silhouette_man)
kmeans_dist_jac<- kmeans(d_matrix_jac, n_clust_silhouette_jac)


similarity_plot(d_matrix_euc, X_sc_kmeans_dist_euc$cluster,X_sc_kmeans_dist_euc$cluster)
similarity_plot(d_matrix_man, X_sc_kmeans_dist_man$cluster,X_sc_kmeans_dist_man$cluster)
similarity_plot(d_matrix_jac, X_sc_kmeans_dist_jac$cluster,X_sc_kmeans_dist_jac$cluster)

#connecting the membership to the player name
est_df<- data.frame(player_slug<- rownames(A_weight), label<-as.numeric(X_sc_kmeans_dist_euc$cluster))
colnames(est_df)<- c('player_slug','est_cl')

# Network plot

#joining the data, inserting covariates
g_tbl<-g %>% as_tbl_graph() %>% activate(edges)%>%select(-.groups)
g_tbl<-g_tbl%>% activate(nodes) %>% inner_join(est_df, by = c('name' = 'player_slug') ) %>% 
  inner_join(players_df, by = c('name' = 'player_slug'))  %>%
  mutate(degree_pl = degree(g_tbl,mode = 'out')/degree(g_tbl,mode = 'all')) %>%
  arrange()


g_tbl <- g_tbl %>%mutate(est_cl = as.factor(est_cl))

g_tbl_df<- g_tbl %>% activate(nodes) %>% as.data.frame()
rank_vs_cluster(g_tbl_df, g_tbl_df$est_cl,est_model = est_model)

#table data---

g_tbl_df%>%group_by(est_cl)%>% summarise(mean(degree_pl))
cor(g_tbl_df$median_rank,as.numeric(g_tbl_df$est_cl))

#----------------------------------------------
#Hierarchical clustering applied to distance matrix

h_clust_euc <-hclust(dist(A_weight,upper = T,method = 'euc'))
h_clust_man <-hclust(dist(A_weight,upper = T,method = 'man'))
h_clust_jac <-hclust(dist(A_weight,upper = T,method = 'bin'))

silhouette_plot_euc<-fviz_nbclust(dist_X_sc_euc, hcut,method = 'sil')
silhouette_plot_man<-fviz_nbclust(dist_X_sc_man, hcut,method = 'sil')
silhouette_plot_jac<-fviz_nbclust(dist_X_sc_jac, hcut,method = 'sil')

n_clust_silhouette_euc<- which.max(silhouette_plot_euc$data$y)
n_clust_silhouette_man<- which.max(silhouette_plot_man$data$y)
n_clust_silhouette_jac<- which.max(silhouette_plot_jac$data$y)

memb_euc<-cutree(h_clust_euc,k = n_clust_silhouette_euc)
memb_man<-cutree(h_clust_man,k = n_clust_silhouette_man)
memb_jac<-cutree(h_clust_jac,k = n_clust_silhouette_jac)


similarity_plot(d_matrix_euc,memb_euc, memb_euc)
similarity_plot(d_matrix_man,memb_man, memb_man)
similarity_plot(d_matrix_jac,memb_jac, memb_jac)

#connecting the membership to the player name
est_df<- data.frame(player_slug<- rownames(A_weight), label<-as.numeric(memb_man))
colnames(est_df)<- c('player_slug','est_cl')

# Network plot

#joining the data, inserting covariates
g_tbl<-g %>% as_tbl_graph() %>% activate(edges)%>%select(-.groups)
g_tbl<-g_tbl%>% activate(nodes) %>% inner_join(est_df, by = c('name' = 'player_slug') ) %>% 
  inner_join(players_df, by = c('name' = 'player_slug'))  %>%
  mutate(degree_pl = degree(g_tbl,mode = 'out')/degree(g_tbl,mode = 'all')) %>%
  arrange()


g_tbl <- g_tbl %>%mutate(est_cl = as.factor(est_cl))


color_palette <- colorRampPalette(c("grey", "black"))(n = 5)
set.seed(123)
ggraph(g_tbl, layout = 'fr') +
  geom_edge_link(aes(alpha = weight)) + # Set custom color palette
  geom_node_point(aes(color = est_cl), size = scale(rowSums(A))+2) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank())
g_tbl_df<- g_tbl %>% activate(nodes) %>% as.data.frame()
rank_vs_cluster(g_tbl_df, g_tbl_df$est_cl,est_model = est_model)

#table data---

g_tbl_df%>%group_by(est_cl)%>% summarise(mean(degree_pl))
cor(g_tbl_df$median_rank,as.numeric(g_tbl_df$est_cl))





n_clust<-which.max(silhouette_vector)+1
print(n_clust)

X_sc_kmeans <- kmeans(X_sc,n_clust)



b<-Spectrum(as.matrix(dist(A,upper = T)))
similarity<- b$eigenvector_analysis



print()

# Random Rewiring
erdos
# Montecarlo approach


fit$memberships
fit$print()
fitted(fit)



