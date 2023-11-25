library(gt)
library(gridExtra)
library(ggplot2)
library(ggside)
library(dplyr)
# monkey application

setwd('/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/monkey application')
monkey_df<-read.csv(file = "dominance.data.csv")
load("df.ranks.Rdata")
df_rank <- get("df.ranks")  

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

df_rank = df_rank %>% filter(name %in% c(e$winner_slug, e$loser_slug))

library(igraph)
library(ggraph)
library(tidygraph)

e_graph = e %>% filter(weight !=0)
g<-graph_from_data_frame(e_graph,directed = T)

g <- g %>% as_tbl_graph()
g%>%activate(nodes) %>% filter(degree(g)<5)

ggraph(g,layout = 'stress')+
  geom_edge_link(aes(alpha = weight), arrow= arrow(length = unit(.25,units = 'cm'),angle = 5))+
  geom_node_point(aes(color = degree(g)))+
  labs(title = 'Network Vervet Monkeys') +
  theme_minimal() 

  
  
  
e = data.frame(e)
colnames(e)<- c('winner_slug','loser_slug','weight')
e= data.frame(e) %>% mutate(winner_slug = as.character(winner_slug)) %>%
  mutate(loser_slug = as.character(loser_slug))
bottom_x <- 1
top_x<- length(unique(c(e$winner_slug, e$loser_slug)))
range_x<- paste0('Range [',bottom_x,'-',top_x,']')

# range_list = list(range1 = c(1,10),
#                   range2= c(21,30),
#                   range3 =c(51,60), 
#                   range4 = c(71,80))
range_list = list(range1 = c(1,length(unique(c(e$winner_slug, e$loser_slug)))))
#storing summary values
summary_df<- data.frame(players=0,probabilities=0,victories=0,games_played=0,range_pl=0)
#storing edgelist
tile_df<- data.frame(player_1=0, player_2=0, probabilities = 0,marginal_x=0,marginal_y=0, range_pl = 0)

for(i in 1:length(range_list)){
  bottom_range<- range_list[[i]][1]
  upper_range<- range_list[[i]][2]
  
  list_pl <- selecting_topx_andbottomy(e,bottom_x,top_x,bottom_range,upper_range)
  
  tile_df_i <- data.frame(player_1= list_pl$e_df$Player_x,
                          player_2= list_pl$e_df$Player_y,
                          probabilities= list_pl$e_df$weight,
                          marginal_x= list_pl$e_df$marginal_x,
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

tile_df1=tile_df

g<-graph_from_data_frame(e,directed = T)

g <- g %>% as_tbl_graph() 
df_rank = g %>% activate(nodes) %>% mutate(rank = degree(g,mode = 'in')/degree(g,mode = 'all')) %>%as.data.frame()


# Heatmap
#-------------------------------------------------------------------------------
tile_df1 = tile_df1 %>% inner_join(df_rank %>% select(name, rank), by = c('player_1' = 'name')) 
tile_df1 = tile_df1 %>% rename(winner_rank = rank)
tile_df1 = tile_df1 %>% left_join(df_rank%>% select(name,rank), by = c('player_2' = 'name') )
tile_df1 = tile_df1 %>% rename(loser_rank = rank)

# Collect the plots in a list
png(paste0("heatmapvs1-90",range_x,".png"),width = 1920, height = 1080)

my_plot<-ggplot(tile_df1, aes(x = reorder(player_2,loser_rank), y = reorder(player_1,winner_rank, decreasing=F))) +
  geom_tile(aes(fill = probabilities), colour = "grey50") +
  scale_fill_gradient(low = "white", high = "red") +
  geom_xsidecol(aes(x = player_2,y = marginal_x, fill = marginal_x))+
  geom_ysidecol(aes(y = player_1,x = marginal_y, fill = marginal_y))+
  labs(title = paste0('Heatmap of ',range_x,' vs Range [1-29]' ), x = "Players in Range [1-15]", y = paste0("Players in Range ", range_x)) +
  theme_bw() +
  theme(legend.direction = "vertical") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

print(my_plot)
