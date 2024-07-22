library(dplyr)
library(ggplot2)

df_matches = read.csv('/Users/lapo_santi/Downloads/archive/atp_matches_till_2022.csv')
df_rank = read.csv('/Users/lapo_santi/Downloads/archive/atp_rankings_till_2022.csv')


# Define the start and end dates for the 2016/2017 season
start_date <- as.Date('2016-12-29')
end_date <- as.Date('2017-11-29')
end_of_season<- as.Date('2017-10-29')

df_rank_sum = df_rank%>%
  mutate(ranking_date = as.Date(as.character(ranking_date), format="%Y%m%d")) %>%
  filter(between(ranking_date,start_date, end_date))%>%
  group_by(player) %>%
  filter(ranking_date == max(ranking_date)) %>%
  select(player,rank)%>%
  rename(last_rank = rank)






matches_2017 = df_matches%>% 
  mutate(tourney_date = as.Date(as.character(tourney_date), format="%Y%m%d"))%>%
  filter(tourney_date >= start_date & tourney_date <= end_date)



name_id_df = matches_2017 %>% select(winner_name, winner_id)%>%
  rename(player = winner_name, id = winner_id)%>%
  rbind(data.frame( player = matches_2017$loser_name, id = matches_2017$loser_id))%>%
  distinct(player,id)




player_list_alphabet = unique(c(matches_2017$winner_name,matches_2017$loser_name))
player_list_alphabet_sorted = data.frame(players = sort(player_list_alphabet), median_rank=NA)

for(p in 1:nrow(player_list_alphabet_sorted)){
  pl_t =  player_list_alphabet_sorted$players[p]
  
  rank_player_t = vector()
  for(t in 1:nrow(matches_2017)){
    
    if(matches_2017$winner_name[t] == pl_t){
      rank_player_t = append(rank_player_t, matches_2017$winner_rank[t])
    }else if(matches_2017$loser_name[t] ==pl_t){
      rank_player_t = append(rank_player_t, matches_2017$loser_rank[t])
    }
    
  }
  player_list_alphabet_sorted$median_rank[p] <- median(rank_player_t)
}

player_list105 = player_list_alphabet_sorted %>% 
  left_join(name_id_df, by=c('players'='player')) %>%
  left_join(df_rank_sum, by = c('id' = 'player')) %>%
filter(last_rank <= 95) 
  

which.max(player_list105$median_rank)

mean_rank = matches_2017 %>%
  filter(winner_name %in% player_list105$players & loser_name %in% player_list105$players)

n= nrow(player_list105)

Y_matrix = matrix(0, nrow = n , ncol = n)

rownames(Y_matrix)<- player_list105$players
colnames(Y_matrix)<- player_list105$players

N_matrix = matrix(0, nrow = n , ncol = n)

rownames(N_matrix)<- player_list105$players
colnames(N_matrix)<- player_list105$players

for(t in 1:nrow(mean_rank)){
  i<- mean_rank$winner_name[t]
  j<- mean_rank$loser_name[t]
  Y_matrix[i,j] = Y_matrix[i,j] + 1
  
  N_matrix[i,j] = N_matrix[i,j] + 1
  N_matrix[j,i] = N_matrix[j,i] + 1
}

which(rowSums(Y_matrix)==0)
which(colSums(Y_matrix)==0)







similarity_plot(Y_matrix, player_list105$median_rank, player_list105$median_rank)


saveRDS(Y_matrix, "./Data/Tennis application/Y_new.rds")
saveRDS(N_matrix, "./Data/Tennis application/N_new.rds")
saveRDS(player_list105,"./Data/Tennis application/median_rank105.rds" )





library(dplyr)
library(ggplot2)
library(tidyverse)


Y_ij <- read.table("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Data/Tennis application/Y_ij.csv",header  = F,row.names = 1,sep = ",")
N_ij <- read.table("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/Data/Tennis application/N_ij.csv",header  = F,row.names = 1,sep = ",")

# Create the data frame
df <- data.frame(players = rownames(Y_ij), win= rowSums(Y_ij)/rowSums(N_ij), loss = colSums(Y_ij)/rowSums(N_ij))

# Gather the data into a tidy format
df_tidy <- df %>%
  tidyr::gather(key = "result", value = "perc", -players)

# Create a separate data frame for ordering
df_order <- df %>%
  dplyr::arrange(desc(wins)) %>%
  dplyr::mutate(order = 1:n())

# Merge the order data frame with the tidy data frame
df_tidy <- df_tidy %>%
  dplyr::left_join(df_order, by = "players")

ggplot(df_tidy, aes(x = reorder(players, order), y = perc, fill = result)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.background = element_rect(fill="white", colour="black", size=0.5,
                                         linetype="solid"),
        legend.key = element_blank()) +
  labs(x = "Players", y = "Percentage", fill = "Result",
       title = "Player Results",
       subtitle = "Ordered by Wins",
       caption = "Source: ATP tournament")
library(dplyr)
library(ggplot2)
library(tidyverse)



# Create the data frame
df <- data.frame(players = rownames(Y_ij), win= rowSums(Y_ij)/rowSums(N_ij))

# Gather the data into a tidy format
df_tidy <- df %>%
  tidyr::gather(key = "result", value = "perc", -players)

# Create a separate data frame for ordering
df_order <- df %>%
  dplyr::arrange(desc(win)) %>%
  dplyr::mutate(order = 1:n())

# Merge the order data frame with the tidy data frame
df_tidy <- df_tidy %>%
  dplyr::left_join(df_order, by = "players")

ggplot(df_tidy, aes(x = reorder(players, order), y = perc)) +
  geom_col(position = "dodge",fill = "green4", color = "black", alpha = 0.7) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.background = element_rect(fill="white", colour="black", size=0.5,
                                         linetype="solid"),
        legend.key = element_blank()) +
  labs(x = "Players", y = "Percentage", fill = "Result",
       title = "Players' win proportions",
       subtitle = "Ordered by win %",
       caption = "Source: ATP tournament")

sd(rowSums(N_ij))

data.frame(games_played = rowSums(N_ij)) %>%
  ggplot(aes(x = games_played)) +
  geom_histogram(binwidth = 5, fill = "green4", color = "black", alpha = 0.7) +
  geom_vline(aes(xintercept = mean(games_played)), color = "red", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = median(games_played)), color = "blue", linetype = "dashed", size = 1) +
  annotate("text", x = mean(data.frame(games_played = rowSums(N_ij))$games_played), 
           y = max(18), label = paste("Mean"), 
           color = "red", angle = 90, vjust = -0.5) +
  annotate("text", x = median(data.frame(games_played = rowSums(N_ij))$games_played), 
           y = max(18), label = paste("Median"), 
           color = "blue", angle = 90, vjust = 1.2) +
  labs(title = "Distribution of Games Played by Players",
       x = "Number of Games Played",
       y = "Frequency") +
  theme_minimal() +
  ylim(0,21)+
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )




