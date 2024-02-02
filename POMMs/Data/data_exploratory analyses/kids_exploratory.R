

setwd('/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/kids application/network/')
e<-read.csv(file = "edges.csv")
e=e[,-4]
e = cbind(e[,2], e[,1], e[,3])
colnames(e) = c('is_liked','likes','weight')
A = matrix(table(e),nrow = 29,ncol = 29)
N= A + lower.tri(A)*t(A) + upper.tri(A)*t(A)

P = A/N
rownames(P) = 1:29
colnames(P) = 1:29

for(i in 1:nrow(e)){
  e$p[i] <- P[cbind(as.numeric(e$winner_slug)[i],as.numeric(e$loser_slug)[i])]
}

for(i in 1:nrow(e)){
  e$marginal_x[i]<- sum(A[e$loser_slug[i],],na.rm = T)/sum(N[e$loser_slug[i],],na.rm = T)
}

for(i in 1:nrow(e)){
  e$marginal_y[i]<- sum(A[,e$winner_slug[i]],na.rm = T)/sum(N[,e$winner_slug[i]],na.rm = T)
}

expand.grid(P)
g<-graph_from_data_frame(e,directed = T)


g <- g %>% as_tbl_graph() %>% activate(nodes) %>% 
  mutate(marginal_y_pl = degree(g,mode = 'out')/degree(g,mode='tot')) %>% 
  mutate(marginal_x_pl = degree(g,mode = 'in')/degree(g,mode='tot'))

pl_df=g%>%activate(nodes)%>%as.data.frame()
new_ord_x = reorder(pl_df$name, pl_df$marginal_x_pl,decreasing = F)
new_ord_y =as.numeric(attributes(new_ord_x)[[1]])

e = e %>% mutate(loser_slug = factor(loser_slug,labels =new_ord_y,ordered = T )) %>% mutate(winner_slug= factor(winner_slug,labels =new_ord_y,ordered = T ))

ggraph(g,layout = 'stress')+
  geom_edge_link(aes(alpha = weight), arrow= arrow(length = unit(.25,units = 'cm'),angle = 5))+
  geom_node_point(aes(color = degree(g)))+
  labs(title = 'Network of kids players', subtitle = 'Who is your best friend?', color = 'number of interactions') +
  theme_minimal() 

ggplot(e, aes(x = loser_slug, y = winner_slug)) +
  geom_tile(aes(fill =p ), colour = "grey50") +
  scale_fill_gradient(low = "white", high = "red") +
  geom_xsidecol(aes(x = winner_slug,y = marginal_x, fill = marginal_x))+
  geom_ysidecol(aes(y = loser_slug,x = marginal_y, fill = marginal_y))+
  labs(title = paste0('Heatmap of ',range_x,' vs Range [1-15]' ), x = "Players in Range [1-15]", y = paste0("Players in Range ", range_x)) +
  theme_bw() +
  theme(legend.direction = "vertical") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
