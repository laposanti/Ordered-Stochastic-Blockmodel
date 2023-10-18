library(gt)

#Uploading data
#-------------------------------------------------------------------------------
match_2017_url <- 'https://pkgstore.datahub.io/sports-data/atp-world-tour-tennis-data/match_scores_2017_unindexed_csv/data/df00561878fee97bf28b92cc70ae1d54/match_scores_2017_unindexed_csv.csv'
ranking_url <- 'https://pkgstore.datahub.io/sports-data/atp-world-tour-tennis-data/rankings_1973-2017_csv/data/79dd58b82401b1e872c23d4d2b6365fb/rankings_1973-2017_csv.csv'

#importing the data

df_rank <- read.csv('/Users/lapo_santi/Desktop/Nial/raw_tennis_data/rankings_1973-2017_csv.csv')

match_2017_url <- 'https://pkgstore.datahub.io/sports-data/atp-world-tour-tennis-data/match_scores_2017_unindexed_csv/data/df00561878fee97bf28b92cc70ae1d54/match_scores_2017_unindexed_csv.csv'
df_match <- read.csv(match_2017_url)

my_edges = df_match %>% select(winner_slug, loser_slug)
my_table <- as.matrix(table(my_edges))
e<-expand.grid(rownames(my_table),colnames(my_table))
for(i in 1:nrow(e)){
  e$weight[i]<- my_table[e$Var1[i], e$Var2[i]]
}

e <- data.frame(winner_slug= e$Var1, loser_slug = e$Var2, weight =e$weight)
e = e %>% inner_join(g_tbl_df %>% select(name, winner_rank), by = c('winner_slug'='name'))
e = e %>% inner_join(g_tbl_df %>% select(name, loser_rank), by = c('loser_slug'='name'))

#input: edgelist
#output: matrix p_ij, dataframe

selecting_topx_andbottomy <- function(e, from_x, to_x, from_y, to_y){
  
  
  
  tot_x <- length(from_x:to_x)
  tot_y <- length(from_y:to_y)
  
  
  from_x_players <- g_tbl_df %>% arrange(median_rank) %>% select(name) %>% head(from_x-1)
  to_x_players <- g_tbl_df %>% arrange(median_rank) %>% select(name) %>% head(to_x)
  from_xto_x_players<- setdiff(unlist(to_x_players), unlist(from_x_players))
  
  from_y_players <- g_tbl_df %>% arrange(median_rank) %>% select(name) %>% head(from_y-1)
  to_y_players <- g_tbl_df %>% arrange(median_rank) %>% select(name) %>% head(to_y)
  from_yto_y_players<- setdiff(unlist(to_y_players), unlist(from_y_players))
  
  A_xvsy <- matrix(0, tot_x, tot_y)
  rownames(A_xvsy) <- unlist(from_xto_x_players)
  colnames(A_xvsy) <- unlist(from_yto_y_players)
  
  for (i in 1:nrow(e)) {
    winner <- e$winner_slug[i]
    loser <- e$loser_slug[i]
    if (winner %in% rownames(A_xvsy) && loser %in% colnames(A_xvsy)) {
      A_xvsy[winner, loser] <- A_xvsy[winner, loser]+ e$weight[i] 
    }
  }
  
  N_xvsy = matrix(0,tot_x,tot_y)
  rownames(N_xvsy) <- unlist(from_xto_x_players)
  colnames(N_xvsy) <- unlist(from_yto_y_players)
  for (i in 1:nrow(e)) {
    winner <- e$winner_slug[i]
    loser <- e$loser_slug[i]
    if (winner %in% rownames(N_xvsy) && loser %in% colnames(N_xvsy)) {
      N_xvsy[winner, loser] <- N_xvsy[winner, loser] + e$weight[i] 
    }
    if(winner %in% colnames(N_xvsy) && loser %in% rownames(N_xvsy)) {
      N_xvsy[loser, winner] <- N_xvsy[loser, winner] + e$weight[i] 
    }
  }
  
  
  p_xvsy <- A_xvsy/N_xvsy
  rownames(p_xvsy) <-  unlist(from_xto_x_players)
  colnames(p_xvsy) <-  unlist(from_yto_y_players)
  
  e <- e %>%
    inner_join(g_tbl_df %>%
                 select(name, median_rank) %>%
                 rename(winner_rank = median_rank),
               by = c('winner_slug' = 'name')) %>%
    inner_join(g_tbl_df %>%
                 select(name, median_rank) %>%
                 rename(loser_rank = median_rank),
               by = c('loser_slug' = 'name'))
  
  e_df<- expand.grid(from_xto_x_players,from_yto_y_players )
  colnames(e_df)<- c('Player_x', 'Player_y')
  
  e_df$weight <- rep(0,nrow(e_df))
  for(i in 1:nrow(e_df)){
    e_df$weight[i]<- p_xvsy[e_df$Player_x[i],e_df$Player_y[i]]
  }
  
  
  e_df$marginal_y<- rep(0,nrow(e_df))
  for(i in 1:nrow(e_df)){
    e_df$marginal_y[i]<- sum(A_xvsy[e_df$Player_x[i],],na.rm = T)/sum(N_xvsy[e_df$Player_x[i],],na.rm = T)
  }
  
  
  
  row_means <- rowSums(A_xvsy,na.rm = T)/(rowSums(N_xvsy,na.rm = T))
  col_means <- colSums(A_xvsy,na.rm = T)/(colSums(N_xvsy,na.rm = T))
  row_df <- data.frame(players_x= names(row_means), row_means = row_means, 
                       col_means=col_means,
                       victories = rowSums(A_xvsy,na.rm = T),
                       games_played =rowSums(N_xvsy,na.rm = T))
  return(list(A_xvsy=A_xvsy, N_xvsy=N_xvsy,e_df= e_df, p_xy = p_xvsy, row_df = row_df))
}

bottom_x <- 1
top_x<- 90
range_x<- paste0('Range [',bottom_x,'-',top_x,']')

# range_list = list(range1 = c(1,10),
#                   range2= c(21,30),
#                   range3 =c(51,60), 
#                   range4 = c(71,80))
range_list = list(range1 = c(1,90))
#storing summary values
summary_df<- data.frame(players=0,probabilities=0,victories=0,games_played=0,range_pl=0)
#storing edgelist
tile_df<- data.frame(player_1=0, player_2=0, probabilities = 0,marginal_y=0, range_pl = 0)

for(i in 1:length(range_list)){
  bottom_range<- range_list[[i]][1]
  upper_range<- range_list[[i]][2]
  
  list_pl <-selecting_topx_andbottomy(e,bottom_x,top_x,bottom_range,upper_range)
  
  tile_df_i <- data.frame(player_1= list_pl$e_df$Player_x,
                          player_2= list_pl$e_df$Player_y,
                          probabilities= list_pl$e_df$weight,
                          marginal_y= list_pl$e_df$marginal_y,
                          range_pl = paste0('Range [',bottom_range,'-',upper_range,']'))
  
  tile_df = rbind(tile_df, tile_df_i)
  
  sum_i <- data.frame(
    players = list_pl$row_df$players_x, 
    probabilities = list_pl$row_df$row_means,
    victories = list_pl$row_df$victories,
    games_played = list_pl$row_df$games_played,
    range_pl = paste0('Range [', bottom_range, '-', upper_range, ']')
  )
  
  summary_df = rbind(summary_df, sum_i)
}
summary_df=summary_df[-1,]
tile_df=tile_df[-1,]

tile_df1=tile_df %>% filter(range_pl == 'Range [1-10]')
tile_df2=tile_df %>% filter(range_pl == 'Range [21-30]')
tile_df3=tile_df %>% filter(range_pl == 'Range [51-60]')
tile_df4=tile_df %>% filter(range_pl == 'Range [71-80]')


# Heatmap
#-------------------------------------------------------------------------------
tile_df1 = tile_df1 %>% inner_join(g_tbl_df %>% select(name, median_rank), by = c('player_1' = 'name')) 
tile_df1 = tile_df1 %>% rename(winner_rank = median_rank)
tile_df1 = tile_df1 %>% left_join(g_tbl_df%>% select(name,median_rank), by = c('player_2' = 'name') )

# Collect the plots in a list
plot_list <- list(
  ggplot(tile_df1, aes(x = reorder(player_2,loser_rank,decreasing = F), y = reorder(player_1,winner_rank,decreasing = T))) +
    geom_tile(aes(fill = probabilities)) +
    scale_fill_gradient(low = "white", high = "red") +
    geom_ysidecol(aes(x = marginal_y, fill = marginal_y)) +
    labs(title = paste0('Heatmap of ',range_x,' vs Range [1-10]' ), x = "Players in Range [1-10]", y = paste0("Players in Range ", range_x)) +
    theme_bw() +
    theme(legend.direction = "vertical") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)),
  
  ggplot(tile_df2, aes(x = player_2, y = player_1)) +
    geom_tile(aes(fill = probabilities)) +
    scale_fill_gradient(low = "white", high = "red") +
    geom_ysidecol(aes(x = marginal_y, fill = marginal_y)) +
    labs(title = paste0('Heatmap of ',range_x,' vs Range [21-30]' ), x = "Players in Range [21-30]", y = paste0("Players in Range ", range_x)) +
    theme_bw() +
    theme(legend.direction = "vertical") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)),
  
  ggplot(tile_df3, aes(x = player_2, y = player_1)) +
    geom_tile(aes(fill = probabilities)) +
    scale_fill_gradient(low = "white", high = "red") +
    geom_ysidecol(aes(x = marginal_y, fill = marginal_y)) +
    labs(title = paste0('Heatmap of ',range_x,' vs Range [51-60]' ), x = "Players in Range [51-60]", y = paste0("Players in Range ", range_x)) +
    theme_bw() +
    theme(legend.direction = "vertical") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)),
  
  ggplot(tile_df4, aes(x = player_2, y = player_1)) +
    geom_tile(aes(fill = probabilities)) +
    scale_fill_gradient(low = "white", high = "red") +
    geom_ysidecol(aes(x = marginal_y, fill = marginal_y)) +
    labs(title = paste0('Heatmap of ',range_x,' vs Range [71-80]' ), x = "Players in Range [71-80]", y = paste0("Players in Range ", range_x)) +
    theme_bw() +
    theme(legend.direction = "vertical") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
)

# Arrange the plots in a 2x2 grid
getwd()
png(paste0("heatmap",range_x,".png"),width = 1920, height = 1080)
grid.arrange(
  plot_list[[1]],
  plot_list[[2]],
  plot_list[[3]],
  plot_list[[4]],
  nrow = 2
)

# Close the device to save the plot
dev.off()
#Bar plot
#-------------------------------------------------------------------------------
png(paste0("barplot",range_x,".png"),width = 817, height = 441)
ggplot(summary_df, aes(x = players, fill = factor(range_pl))) +
  geom_col(aes(y = factor(games_played), fill = 'Games Played'), position = "identity") +
  geom_col(aes(y = factor(victories),fill = 'Victories'), position = "identity") +
  facet_wrap(~ factor(range_pl), nrow = 1) +
  labs(title = paste0("Victories and Games Played of Players", range_x," vs. Different Ranges"),
       x = "Players",
       y = "Count",
       fill = "Data Type") +
  scale_fill_discrete(name = "Data Type",
                      labels = c("Games Played", "Victories")) +
  ylim(factor(0:max(summary_df$games_played)))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()


#Box plot
#-------------------------------------------------------------------------------
png(paste0("boxplot",range_x,".png"),width = 817, height = 441)
ggplot(summary_df, aes( y = probabilities, fill = range_pl)) +
  geom_boxplot() +
  facet_wrap(~ factor(range_pl),nrow = 1) +
  labs(title = paste0("Victories of Players", range_x," vs. Different Ranges"),
       x = "Players",
       y = "Probability of Victory",
       fill = "Player Range") +
  theme_minimal()
dev.off()
#Table
#-------------------------------------------------------------------------------


sum_table<-summary_df%>%group_by(range_pl)%>% summarise(mean_games = median(games_played,na.rm=T),
                                                        mean_victories = median(victories,na.rm=T),
                                                        median = median(probabilities,na.rm=T),
                                                        inter= IQR(probabilities,na.rm=T),
                                                        sd = sd(probabilities,na.rm=T))

setwd('/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/exploratory_graphs_with_inference/exploratoryII/ExploratoryIII')
save_table_to_file(sum_table,filename =paste0(range_x,'.csv'),title = paste0(range_x,'.csv') )




