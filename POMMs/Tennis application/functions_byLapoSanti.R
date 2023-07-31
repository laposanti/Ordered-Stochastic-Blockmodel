####################################################################################
# SPECTRAL CLUSTERING ####
####################################################################################
#Input:
#adj <-  the adjacency matrix of the graph
#K is the number of clusters
spectral_clustering <- function(adj, K){
  #computing the degree vector
  d <- rowSums(adj)
  #computing the laplacian
  lapl <- diag(d) - adj
  #eigenvalues of the laplacian
  ev <- eigen(lapl, symmetric = TRUE)
  #extracting K  eigenvectors
  K_eigenvectors <- ev$vectors[ , (ncol(ev$vectors)-K+1):ncol(ev$vectors)] 
  #applying the kmeans algorithm
  ee <- kmeans(K_eigenvectors,K)$cluster
  return(ee)
}

####################################################################################
# Graph PLOT ####
####################################################################################
#z_vector is a membership vector of lenth n
#degree_vector is a vector containing as entries each node's degree
plotting_graphs <- function(adj, z_vector, degree_vector){
  g <- graph_from_adjacency_matrix(adj)
  ggraph(g,layout = "fr")+
    geom_edge_link0(edge_colour = "gray")+
    geom_node_point(aes(fill =as.character(z_vector ), size=degree_vector ),shape = 21)+
    scale_edge_width(range = c(0.2,3))+
    scale_colour_gradient2(low = "red", mid = "white", high = "blue") +
    scale_size(range = c(1,6))+
    theme_graph()+
    theme(legend.position = "none")
}

plotting_centrality <- function(a, nodes_labels, centrality_measure, custom_title, perclabels){
  labels_prva <- nodes_labels
  p1 <- ggraph(a,layout = "centrality", centrality = centrality_measure, tseq = seq(0,1,0.15)) +
    draw_circle(use = "cent") +
    geom_edge_link0(edge_color="gray",edge_width=0.3)+
    geom_node_point(aes(fill=as.factor(z_0)),size=2,shape=21)+
    scale_colour_gradient(low = "red", high = "blue")+
    geom_node_text(aes(filter = centrality >= quantile(centrality_measure, probs = (perclabels)), label = labels_prva), family = "serif") +
    theme_graph()+
    theme(legend.position = "none")+
    coord_fixed()+
    labs(title= custom_title)
  print(p1)
}

####################################################################################
# from classical Edgelist to 3 col matrix with 0-1 on the third col  ####
####################################################################################
#

proper_matrix_form <- function(edgelist){
  if(is.data.frame(edgelist)){
    colnames(edgelist) <- c("x1","x2")
    row.names(x) <- c(1:dim(edgelist)[1])
    #all nodes id together
    v_nodes <- append(v_nodes,edgelist$x1)
    v_nodes <- append(v_nodes,edgelist$x2)
    #number of nodes
    num_nodes <- as.numeric(length(unique(v_nodes)))
    #unique id of nodes
    v_unique <- unique(v_nodes)
    #creating 3 support vector
    vector_support <- matrix(NA, nrow = dim(edgelist)[1], ncol= 2)
    for(i in 1:dim(edgelist)[1]){
      vector_support[i,1] <- edgelist[[i,1]]
      vector_support[i,2] <- edgelist[[i,2]]
    }
    
    #vettore contente 1,2,3,...,516 x 516 volte
    
    vector_support_1 <- rep.int(v_unique, num_nodes)
    
    #vettore contente 1 x 516 volte, 2 x 516 volte, ..., 516 x 516 volte
    
    vector_support_2 <- vector()
    for(i in 1:num_nodes){
      vector_support_2 <- append(vector_support_2, rep.int(v_unique[i],num_nodes))
    }
    
    #matrice contenente le 516^2 possibili combinazioni
    
    matrix_1 <- matrix(NA, nrow = length(vector_support_2), ncol=3 )
    for(i in 1:length(vector_support_2)){
      matrix_1[i,1] <- vector_support_1[i]
      matrix_1[i,2] <- vector_support_2[i]
    }
    for(i in 1:nrow(matrix_1)){
      j=1
      find_link = 0
      while(j<=nrow(vector_support) && find_link == 0){
        if(matrix_1[i,1] == vector_support[j,1] && matrix_1[i,2] == vector_support[j,2]){
          matrix_1[i,3] <-  1
          find_link = 1
        }
        else if(matrix_1[i,1] == vector_support[j,2] && matrix_1[i,2] == vector_support[j,1]){
          matrix_1[i,3] <-  1
          find_link = 1
        }
        j=j+1
      }
      if(find_link == 0){
        matrix_1[i,3] <- 0
      }
      if (i%%10000 == 0){print(paste("Iteration:", i))}
    }
  }
  else{
    print("x must be a dataframe")
  }
  print(matrix_1)
  return(matrix_1)
}


####################################################################################
# ADJACENCY MATRIX PLOT ####
####################################################################################
#Input:
# adj <-adjacency matrix of the network
# z_0 <- any membership vector

adjacency_plot <- function(adj, z_0){
  ID <- c(1:nrow(adj))
  df_a <- data.frame(adj)
  #sorting columns
  z_0_a <- data.frame(z_0, ID )
  df_a <- data.frame(adj, ID)
  df_a <- left_join(df_a, z_0_a, by="ID")
  df_a <- df_a[order(df_a$z_0, decreasing = FALSE),]
  df_a$z_0 <- NULL
  df_a$ID <- NULL
  #sorting rows
  df_a <- t(df_a)
  df_a <- data.frame(df_a, ID)
  df_a <- left_join(df_a, z_0_a, by="ID")
  df_a <- df_a[order(df_a$z_0),]
  df_a$z_0 <- NULL
  df_a$ID <- NULL
  
  a_aux <- as.matrix(df_a)
  h_aux <- heatmap(a_aux , col=0:1 , symm = TRUE, Rowv = NA, Colv = NA, ColSideColors = as.character(z_0[order(z_0)]), RowSideColors = as.character(z_0[order(z_0, decreasing = FALSE)]) )
  return(h_aux)
}
####################################################################################
# figure 3 ploT ####
####################################################################################


adj_durante <- function(Y,z_0){
  a_prova <- round(length(unique(z_0))/3)
  b_prova <- round(length(unique(z_0))/3)
  c_prova <- length(unique(z_0)) - a_prova - b_prova
  
  Y_sup <- Y
  diag(Y_sup) <- 0
  V <- nrow(Y_sup)
  row_plot_Y_sup <- as.data.frame(as.factor(matrix((z_0),V,1)))
  names(row_plot_Y_sup) <-"z_0"
  rownames(Y_sup) <- rownames(row_plot_Y_sup)
  mycolors <- c(brewer.pal(10,"RdBu")[c(1:a_prova)],brewer.pal(10,"PRGn")[c(1:b_prova)],brewer.pal(9,"YlOrBr")[c(1:c_prova)])
  names(mycolors) <- unique(row_plot_Y_sup$z_0)
  mycolors <- list(z_0 = mycolors)
  
  Network <- pheatmap(Y_sup,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_Y_sup,annotation_names_row=F,show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)
  
  # ------------------------------------
  
  g <- grid.arrange(Network[[4]],nrow=1,vp=viewport(width=1, height=1))
  g2 <- cowplot::ggdraw(g)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))
  
  print(g2)
}

####################################################################################
# figure 4 ploT ####
####################################################################################

figure_4 <- function(sim_matrix, z_clustering){
V <- length(z_clustering)
sim_matrix <-   sm_fs
z_clustering <- z_vi_fs
a_prova <- round(length(unique(z_clustering))/3)
b_prova <- round(length(unique(z_clustering))/3)
c_prova <- length(unique(z_clustering)) - a_prova - b_prova

row_plot_GN <- as.data.frame(as.factor(matrix(z_clustering,V,1)))
names(row_plot_GN) <- "z_clustering"
rownames(sim_matrix) <- rownames(row_plot_GN)
mycolors <- c(brewer.pal(10,"RdBu")[c(1:a_prova)],brewer.pal(10,"PRGn")[c(1:b_prova)],brewer.pal(9,"YlOrBr")[c(1:c_prova)])
names(mycolors) <- unique(row_plot_GN$z_clustering)
mycolors <- list(z_clustering = mycolors)

Marg <- pheatmap(sim_matrix,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_GN,annotation_names_row=F,show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)

# ------------------------------------

sim_matrix_x <- pr_cc(Z_GN_x[,(burn_in+1):N_iter])
z_clustering_VI_x <- minVI(sim_matrix_x,method="avg",max.k=20)
z_clustering_x <- z_clustering_VI_x$cl

row_plot_GN <- as.data.frame(as.factor(matrix(z_clustering_x,V,1)))
names(row_plot_GN) <- "z_clustering_x"
rownames(sim_matrix_x) <- rownames(row_plot_GN)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3])
names(mycolors) <- unique(row_plot_GN$z_clustering_x)
mycolors <- list(z_clustering_x = mycolors)

Cov <- pheatmap(sim_matrix_x,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_GN,annotation_names_row=F,show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)

# ------------------------------------

g <- grid.arrange(Marg[[4]],Cov[[4]],ncol=1,vp=viewport(width=1, height=1))
g2 <- cowplot::ggdraw(g)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))

print(g2)
}

####################################################################################
# Clusters' HH index PLOT ####
####################################################################################

#computing concentration index for each cluster

plotting_clusters_HHindex <- function(Y, z_est){
  z_est <- z_est
  K <- nrow(unique(z_est))
  
  HH_matrix <- matrix(NA, nrow= K , ncol=2)
  
  for(i in unique(z_est)){
    classes <-  y$Disease.class[which(z_est==i)]
    print(length(classes))
    counts <- table(classes)
    HH <- 0
    for (j in 1:length(counts)){
      HH <- HH + (counts[j]/sum(counts))^2
    }
    HH_matrix[i,1] <- as.numeric(HH)
    HH_matrix[i,2] <- sum(counts)/length(unique(counts))*0.01
  }
  HH_matrix <- as.data.frame(HH_matrix)
  HH_matrix <- HH_matrix[order(as.numeric(rownames(HH_matrix)), decreasing = FALSE),]
  
  ggplot(HH_matrix, aes(x=as.numeric(rownames(HH_matrix)) , y = HH_matrix[,c(1,2)])) + 
    geom_line(aes(y=HH_matrix[,1], color="Actual concentration")) + 
    geom_line(aes(y=HH_matrix[,2], color="Equal concentration")) +
    labs(y = "Concentration", x = "clusters")+
    scale_x_discrete(limits=as.character(as.numeric(rownames(HH_matrix))))
  
}


####################################################################################
# Clusters' degree PLOT ####
####################################################################################

plotting_clusters_degree <- function(Y, z_est){
  n <- dim(Y)[1]
  g <- graph_from_adjacency_matrix(Y)
  m <- dim(as_edgelist(g))[1]/2
  K <- nrow(unique(z_est))
  
  erdos <- erdos.renyi.game(n, m, type="gnm") 
  erdos_adj <- as_adjacency_matrix(erdos)
  nodes_degrees_erdos <- rowSums(erdos_adj) 
  
  
  Degree_matrix <- matrix(NA, nrow= K, ncol=2)
  nodes_degrees <- rowSums(Y)
  
  for(i in unique(z_est)){
    degrees <-  nodes_degrees[which(z_est==i)]
    total_degree <- sum(degrees)/sum(z_est==i)
    
    erdos_degree <- nodes_degrees_erdos[which(z_est==i)]
    total_erdos_degree <- sum(erdos_degree)/sum(z_est==i)
    
    Degree_matrix[i,1] <- total_degree
    Degree_matrix[i,2] <- total_erdos_degree
  }
  
  Degree_matrix <- as.data.frame(Degree_matrix)
  Degree_matrix <- Degree_matrix[order(as.numeric(rownames(Degree_matrix)), decreasing = FALSE),]
  
  my_graph <- ggplot(Degree_matrix, aes(x=as.numeric(rownames(Degree_matrix)) , y = Degree_matrix[,c(1,2)])) + 
    geom_line(aes(y=Degree_matrix[,1], color="Actual cluster degree")) + 
    geom_line(aes(y=Degree_matrix[,2], color="Equal cluster degree")) +
    labs(y = "Degree", x = "clusters")+
    scale_x_discrete(limits=as.character(as.numeric(rownames(Degree_matrix))))
  return(my_graph)
}








####################################################################################
# SIMILARITY MATRIX PLOT ####
####################################################################################
#Input:
# adj <-adjacency matrix of the network/ similarity matrix
# z_0 < true clustering
# z_est <- estimated partition 
similarity_plot <- function(adj, z_0, z_est){
  ID <- c(1:dim(adj)[1])
  df_a <- data.frame(adj)
  #sorting columns accordin to z_0 
  #(descending order to have the clusters along the main diagonal)
  z_0_a <- data.frame(z_0, ID )
  df_a <- data.frame(adj, ID)
  df_a <- left_join(df_a, z_0_a, by="ID")
  df_a <- df_a[order(df_a$z_0, decreasing = FALSE),]
  df_a$z_0 <- NULL
  df_a$ID <- NULL
  #sorting rows according to z_0
  df_a <- t(df_a)
  df_a <- data.frame(df_a, ID)
  df_a <- left_join(df_a, z_0_a, by="ID")
  df_a <- df_a[order(df_a$z_0),]
  df_a$z_0 <- NULL
  df_a$ID <- NULL
  
  #ordering the nodes of z_est for consistency with z_0
  z_aux <- (z_0_a[order(z_0),])
  df_est_aux <- data.frame(z_est, ID)
  # to plot the membership on the columns
  z_est_ord <- vector()
  for(i in 1:nrow(z_aux)){
    z_est_ord <- append(z_est_ord, z_est[z_aux$ID[i]])
  }
  
  
  #to plot the membership on the rows (descending order)
  z_est_desc <- vector()
  for(i in 1:nrow(z_aux)){
    z_est_desc <- append(z_est_desc, z_est[z_aux$ID[i]])
  }
  
  
  a_aux <- as.matrix(df_a)
  h_aux <- heatmap(a_aux , col=0:1 , symm = TRUE, Rowv = NA, Colv = NA, ColSideColors = as.character(z_est_ord), RowSideColors = as.character(z_est_desc))
  return(h_aux)
}
#######################################################################
# Mixed_membership PLOT ####
#######################################################################
mixedmembership_plot <- function(adj, z_0, z_est){
  ID <- c(1:dim(adj)[1])
  df_a <- data.frame(adj)
  #sorting columns accordin to z_0 
  #(descending order to have the clusters along the main diagonal)
  z_0_a <- data.frame(z_0, ID )
  df_a <- data.frame(adj, ID)
  df_a <- left_join(df_a, z_0_a, by="ID")
  df_a <- df_a[order(df_a$z_0, decreasing = FALSE),]
  df_a$z_0 <- NULL
  df_a$ID <- NULL
  #sorting rows according to z_0
  df_a <- t(df_a)
  df_a <- data.frame(df_a, ID)
  df_a <- left_join(df_a, z_0_a, by="ID")
  df_a <- df_a[order(df_a$z_0),]
  df_a$z_0 <- NULL
  df_a$ID <- NULL
  
  #ordering the nodes of z_est for consistency with z_0
  z_aux <- (z_0_a[order(z_0),])
  df_est_aux <- data.frame(z_est, ID)
  # to plot the membership on the columns
  z_est_ord <- vector()
  for(i in 1:nrow(z_aux)){
    z_est_ord <- append(z_est_ord, z_est[z_aux$ID[i]])
  }
  
  
  #to plot the membership on the rows (descending order)
  z_est_desc <- vector()
  for(i in 1:nrow(z_aux)){
    z_est_desc <- append(z_est_desc, z_est[z_aux$ID[i]])
  }
  
  
  a_aux <- as.matrix(df_a)
  h_aux <- heatmap(a_aux , symm = TRUE, Rowv = NA, Colv = NA, ColSideColors = as.character(z_est_ord), RowSideColors = as.character(z_est_desc))
  return(h_aux)
}
  #######################################################################
  # TRACE PLOT ####
  #######################################################################
  #Diagnostic tool to check MCMC mixing behavior
  #it needs a matrix where the columns are the parameters
  # and the rows are the simulated samples
  
  traceplot_lapo <- function(MCMC){
    if("bayesplot" %in% (.packages())){
      a_mcmc <- t(MCMC)
      colnames(a_mcmc) <- c(1:nrow(MCMC))
      mcmc_trace(a_mcmc, pars= c("1","2","3","4","5", "6"))
    }
    else{
      print("Install bayesplot package")
    }
  }
  
  
  
####################################################################################
# APPROXIMATING MARGINAL log-LIKELIHOOD VIA HARMONIC MEAN ####
####################################################################################
#MCMC_sample is the output of the GIbbs Sampler

#returns the negative marginal log-likelihood
marginal_likelihood <- function(adj,a,b,MCMC_sample){
  empty <- vector()
  for(i in 1:ncol(MCMC_sample)){
    if (i%%10000 == 0){print(paste("Iteration:", t))}
    log_likelihood <- log_pY_z(adj,MCMC_sample[, i],a, b)
    likelihood <- 1/exp(log_likelihood)
    empty <- append(empty,likelihood )
  }
  marginal <- 1/(sum(empty)/ncol(MCMC_sample))
  log_marginal <- log(marginal)
  return(log_marginal)
}

marginal_likelihood_1 <- function(adj,a,b,MCMC_sample){
  
  empty <- matrix(NA, nrow=ncol(MCMC_sample), ncol=1)
  
  for(t in 1:ncol(MCMC_sample)){
    log_likelihood <- log_pY_z(adj,MCMC_sample[, t],1, 1)
    likelihood <- 1/exp(log_likelihood)
    empty[t] <- likelihood
    if (t%%10000 == 0){print(paste("Iteration:", t))}
  }
  
  marginal <- 1/(sum(empty)/ncol(MCMC_sample))
  log_marginal <- log(marginal)
}

####################################################################################
# EXPECTED LOSS UNDER VI WHEN THE TRUE PARTITION IS KNOWN ####
####################################################################################
#z0 is the true partition
#MCMC_sample is the output of the GIbbs Sampler

expected_loss_VI <- function(MCMC_sample, z0){
  if (!is.vector(z0)){
    stop("The provided true partition is not a vector")
  }
  else{
    empty_vector <- vector()
    #computing VI between z0 and each MCMC posterior partition
    
    for(i in 1:ncol(MCMC_sample)){
      loss_IV <- vi.dist(z0,MCMC_sample[,i])
      empty_vector <- append(empty_vector,loss_IV)
    }
    
    empty_vector <- as.vector(empty_vector)
    #computing the expected value
    sum_sym <- t(empty_vector)%*%matrix(1,nrow = ncol(MCMC_sample),ncol= 1)
    n_sym <- t(matrix(1,nrow = ncol(t(empty_vector)),ncol= 1))%*%matrix(1,nrow = ncol(t(empty_vector)),ncol= 1)
    posterior_average <- as.numeric(sum_sym)/as.numeric(n_sym)
    return(posterior_average)
  }
}
####################################################################################
# ESTIMATING THETA ####
####################################################################################

thetamatrix <- function(memb,Y,a,b){
  # in: vector of cluster labels (memb), VxV adjancency matrix (Y) and hyperparameters beta priors (a,b)
  # out: vector of Bernoulli log-likelihoods for the edges under ESBM (conditioned on memb and on block-probabilities)

z <- dummy(memb)
H <- ncol(z)
Abs_Freq <- t(z)%*%Y%*%z
diag(Abs_Freq) <- diag(Abs_Freq)/2
Tot <- t(z)%*%matrix(1,V,V)%*%z
diag(Tot) <- (diag(Tot)-table(memb))/2
Block_freq <- (a+Abs_Freq)/(a+b+Tot)
return(Block_freq)
}



####################################################################################
# ESTIMATING THE NUMBER OF CLUSTERS AND QUANTILES ####
####################################################################################
#this function returns the posterior number of clusters
#along with the 1-alfa/2, a/2 quantiles

  K_posterior <- function(MCMC_sample, alfa){
    aux <- vector()
    N_iter <- ncol(MCMC_sample)
    for(i in 1:N_iter){
      z <- as.vector(unique(MCMC_sample[,i]))
      number_communities <- length(z)
      aux <- append(aux,number_communities)
    }
    estimated_number <- sum(aux)/N_iter
    upperquant <- quantile(aux, probs = (1-(alfa/2)))
    lowerquant <- quantile(aux, probs = (alfa/2))
    histogram <- aux
    c <- list(estimated_number, upperquant, lowerquant, histogram)
    return(c)
  }
########################################################
#MIXED MEMBERSHIP STOCHASTIC BLOCK MODELS
########################################################
#-------------------------------------------------------
#                GENERATING A NETWORK
#-------------------------------------------------------

#DEPRECATED

# simulating_mmsbm <- function(N,K,my_eta, my_alpha){
#   #this matrix contains a number N of k-dimensional mixed membership vectors pi_p
#   pi <- matrix(NA, nrow = K, ncol=N)
#   for(i in 1:N){
#     #sampling from a dirichlet random vectors for each node p
#     pi[,] <-  rdirichlet(n = 1, alpha = my_alpha)
#   }
#   #this matrix contains the membership of vector p when interacting with vector q
#   membership_matrix <- matrix(NA, nrow = N, ncol=N)  
#   #this is the final adjacency matrix
#   my_graph <- matrix(NA, nrow=N, ncol=N)
#   for(i in 1:N){
#     for (j in 1:i-1){
#       if(j == 0){
#         my_graph[i, j] <-  0
#       }
#       else{
#         #sampling from the categorical with the pre-specified probabilities
#         z_ij <- rmultinom(n = 1,size = 1, prob = pi[,i])
#         z_ji <- rmultinom(n = 1,size = 1, prob = pi[,j])
#         membership_matrix[i,j] <- which(z_ij >0)
#         membership_matrix[j,i] <- which(z_ji >0)
#         #sampling from a bernoulli of parameter z_ij*eta*z_ji
#         my_graph[i, j] <-  rbinom(1,1,my_eta[which(z_ij >0),which(z_ji >0)])
#         #enforcing simmetry in the adjacency matrix
#         my_graph[upper.tri(my_graph)] = t(my_graph)[upper.tri(my_graph)]
#       }
#     }
#   }
#   
#   #enforcing simmetry in the membership vector
#   membership_matrix[upper.tri(membership_matrix)] = t(membership_matrix)[upper.tri(membership_matrix)]
#   diag(my_graph) <- 0
#   return(my_graph)
# }
#-------------------------------------------------------
#                GENERATING A NETWORK pt-2
#-------------------------------------------------------


simulating_mmsbm <- function(N,K,C, my_alpha, directed){
  C <- as.matrix(s)
  #this matrix contains a number N of k-dimensional mixed membership vectors pi_p
  pi <- matrix(NA, nrow = K, ncol=N)
  for(i in 1:N){
    #sampling from a dirichlet random vectors for each node p
    pi[,i] <-  sample_dirichlet(n=1, alpha = myalpha)
  }
  
  #this matrix contains the membership of vector p 
  membership_matrix <- matrix(NA, nrow = K, ncol=1)  
  
  #this is the final adjacency matrix
  my_graph <- matrix(NA, nrow=N, ncol=N)
  
  #getting the nodes'clusters matrix indicator
  for(i in 1:N){  
    for(j in 1:N){
      #sampling from the categorical with the pre-specified probabilities
      z_ij <- rmultinom(n = 1,size = 1, prob = pi[,i])
      z_ij <- as.matrix(z_ij)
      z_ij_1 <- rmultinom(n = 1,size = 1, prob = pi[,j])
      #storing the results to knwon the true membership vector
      membership_matrix[i] <- which(z_ij>0)
      #sampling from a bernoulli of parameter p <- z_ij%*%C%*%t(z_ij)
      my_graph[i, j] <-  rbinom(1,1,t(z_ij)%*%C%*%z_ij_1)
      }
    
  }
  if(directed == FALSE){
    my_graph[lower.tri(my_graph)] = t(my_graph)[lower.tri(my_graph)]
  }
  results <- list("my_MMSBM_graph" = my_graph, "true membership" = membership_matrix)
  return(results)
}

#-------------------------------------------------------
#                GIBBS SAMPLER
#-------------------------------------------------------

mmsbm.gibbs.lapo <- function(Y, K, N_iter, my_alfa, a,b){

  
  #initialize C
  my_C<- matrix(NA, nrow=K, ncol= K)
  for(i in 1:K){
    for(j in 1:K){
      my_C[i,j] <- rbeta(1, a, b)
    }
  }
  
  #initialize pi
  pi <- matrix(NA, nrow = K, ncol=N)
  for(i in 1:N){
    #sampling from a dirichlet random vectors for each node p
    pi[,] <-  rdirichlet(n = 1, alpha = my_alpha)
  }
  #initialize S_ij
  
  S_ij <- array(NA, c(5,N,N))
  for(i in 1:N){
    for(j in 1:N){
      S_ij[,i,j] <- rmultinom(n = 1,size = 1, prob = pi[,j])
    }
  }
  
  R_ij <- array(NA, c(5,N,N))
  for(i in 1:N){
    for(j in 1:N){
      R_ij[,i,j] <- rmultinom(n = 1,size = 1, prob = pi[,j])
    }
  }
  
  
  
  #creating the empty containers for the samples
  K*K
  C_MCMC <- array(NA, c(K,K,N_iter))
  pi_MCMC <- matrix(NA, K, ncol=N_iter)
  Sij_MCMC <- array(NA, c(K,N,N,N_iter))
  Rij_MCMC <- array(NA, c(K,N,N,N_iter))
  
  #sampling procedure
  for(t in 1:N_iter){
    for(i in 1:N){
      #updating pi
      pi_MCMC[,t]<- my_alpha + rowSums(S_ij, dims=2)[,i] +rowSums(R_ij, dims=2)[,i] 
      
      #updating Sij
      for(j in 1:N){
        p <- matrix(NA, nrow=K)
        for(k in 1:K){
          p[k] <- pi_MCMC[k]*(my_C[k,which(R_ij[,i,j]>0)]^Y[i,j])*(1-my_C[k,which(R_ij[,i,j]>0)])^(1-Y[i,j])
        }
        Sij <- rmultinom(1,1,p)
        Sij_MCMC[,i,j,t] <- Sij
      }
      
      #updating R_ij
      for(j in 1:N){
        q <- matrix(NA, nrow=K)
        for(k in 1:K){
          q[k] <- pi_MCMC[k]*(my_C[k,which(S_ij[,i,j]>0)]^Y[i,j])*(1-my_C[k,which(S_ij[,i,j]>0)])^(1-Y[i,j])
        }
        Rij <- rmultinom(1,1,p)
        Rij_MCMC[,i,j,t] <- Rij
      }
      #updating C (beta generated matrix)
      
      empty_p <- vector()
      empty_q <- vector()
      for(p in 1:N){
        for (q in 1:N){
          if (sum(c(S_ij[,p,q] == S_ij[,1,1]))==K && sum(c(R_ij[,p,q] == R_ij[,1,1]))==K){
            empty_p <- append(empty_p,p)
            empty_q <- append(empty_q,q)
          }
          
        }
      }
      
      #counting number of links among clusters i,j
      sum_y_ij <- 0 
      for(i in 1:length(empty_p)){
        if(my_graph[empty_p[i],empty_q[i]] ==1){
          sum_y_ij <- sum_y_ij + 1
        }
      }
      
      a_l <- my_C[which(S_ij[,1,1]>0)] + sum_y_ij
      b_l <- my_C[which(R_ij[,1,1]>0)]+ sum_y_ij
      my_C[which(S_ij[,1,1]>0),which(R_ij[,1,1]>0)] <- rbeta(1,a_l,b_l)
      C_MCMC[,,t] <- my_C
    }
  }
  return(Sij_MCMC)
}

#-----------------------------------------------------------
#         SUBSETTING GRAPH INDICES FOR CROSS VALIDATION
#-----------------------------------------------------------
# my_data <- adjacency matrix
# K <- number of splits in the dataset
# my_perc <- percentage of data composing the test set

# This functions returns two matrices. The first one contains the nodes'
# positions for the test set and the second one for the training set.
# Each row for the two matrices are sets of indices:
# 1) Train the data on the nodes indicated by the positions of the i-th row of the second matrix;
# 2) Test the data on the nodes indicated by the positions of the i-th row of the first matrix.


creating_cv_sets <- function(my_data, K){
  n <- dim(my_data)[1]
  my_perc = K/n
  
  subsetted_list <- c(1:n)
  subsetted_list_train <- c(1:n)
  
  test_id_matrix <- matrix(NA, nrow=K, ncol = round(my_perc*n, digits=0))
  train_id_matrix <- matrix(NA, nrow=K,ncol= n - round(my_perc*n, digits=0))
  
  for(i in 1:K){
    test_id <- round(runif(length(subsetted_list), min = 0, max=100), digits=3)
    test_id <- as.vector(test_id)
    my_tail <- as.vector(tail(sort(test_id), my_perc*n))
    indices <- vector()
    for(j in 1:length(test_id)){
      for(k in 1:length(my_tail)){
        if(test_id[j] == my_tail[k])
          indices <-  append(indices, j)
      }
    }
    
    aux <- matrix(NA, nrow=length(subsetted_list), ncol=1)
    for(p in 1:length(indices)){
      aux[indices[p]] <- 1
    }
    aux[is.na(aux)] <- 0
    
    final_index <- diag(as.list(aux), nrow=length(subsetted_list)) %*% subsetted_list
    final_index <- subset(final_index, final_index>0)
    test_id_matrix[i,] <- final_index
    train_id_matrix[i,] <- subsetted_list_train[-final_index]
    subsetted_list <- subsetted_list[-indices]
    print(length(subsetted_list))
  }
  result <- list("test_set" = test_id_matrix,"train_set" =  train_id_matrix)
  return(result)
}








#-------------------------------------------------------
#                GIBBS SAMPLER pt 2
#-------------------------------------------------------



mmsbm.gibbs.lapo <- function(Y, K, N_iter, my_alpha, a,b){
  
  
  #ordinare i dati in maniera posizionale
  Y
  N <- dim(Y)[1]
  K <- 5
  a <- 1
  b <- 1
  #initialize C
  my_C<- matrix(NA, nrow=K, ncol= K)
  for(i in 1:K){
    for(j in 1:K){
      my_C[i,j] <- rbeta(1, a, b)
    }
  }
  i <- 1
  j <- 1
  
  #initialize pi and Z_ij
  
  pi <- matrix(NA, nrow = K, ncol=N)
  Z_ij <- matrix(NA, nrow= K,ncol= N)
  for(i in 1:N){
    #sampling from a dirichlet random vectors for each node p
    pi[,] <-  rdirichlet(n = 1, alpha = my_alpha)
    
    for(j in 1:N){
      S_ij[,i,j] <- rmultinom(n = 1,size = 1, prob = pi[,j])
    }
  }
  i <- 1
  j <- 1
  
  R_ij <- array(NA, c(5,N,N))
  for(i in 1:N){
    for(j in 1:N){
      R_ij[,i,j] <- rmultinom(n = 1,size = 1, prob = pi[,j])
    }
  }
  
  
  #creating the empty containers for the samples
  
  C_MCMC <- array(NA, c(K,K,N_iter))
  pi_MCMC <- matrix(NA, K, ncol=N_iter)
  Sij_MCMC <- array(NA, c(K,N,N,N_iter))
  # Rij_MCMC <- array(NA, c(K,N,N,N_iter))
  
  #------------------------------------
  #Beginning of the Iterative procedure
  #------------------------------------
  for(t in 1:N_iter){
    for(i in 1:N){
      
      #updating pi
      pi <- my_alpha + rowSums(S_ij, dims=2)[,i] +rowSums(S_ij, dims=2)[,i]
      pi_MCMC[,t]<- pi
      
      #considering the stochastic update for each node j
      for(j in 1:N){
        
        #defining the multinomial probability
        p <- matrix(NA, nrow=K, ncol=1)
        #sampling the probability for each cluster
        for(k in 1:K){
          p[k] <- pi[k]*(my_C[k,which(S_ij[,i,j]>0)]^Y[i,j])*(1-my_C[k,which(R_ij[,i,j]>0)])^(1-Y[i,j])
        }
        
        #sampling the membership of node i when interacting with node j
        #memorizza l'id del cluster
        S_ij[,i,j] <- which(rmultinom(1,1,p)>0)
        Sij_MCMC[,i,j,t] <- S_ij[,i,j]
        
        #cambiare notazione da Rij a Rji
        #updating Rij as the transpose of Sij
        for(h in 1:K){
          R_ij[h,,] <- t(S_ij[h,,])
        }
        aux_1 <- R_ij[,i,j]
        
        # updating C (beta generated matrix)
        empty_p <- vector()
        empty_q <- vector()
        
        # selecting the indeces of nodes belonging to block K
        
        for(el_1 in 1:N){
          for (el_2 in 1:N){
            if (sum(c(S_ij[,el_1,el_2] == S_ij[,i,j]))==K && sum(c(R_ij[,el_1,el_2] == R_ij[,i,j]))==K){
              empty_p <- append(empty_p,el_1)
              empty_q <- append(empty_q,el_2)
            }
          }
        }
        
        #selecting the indeces
        
        sum_y_ij <- 0 
        for(index_i in 1:N){
          for(index_j in 1:N){
            if(new_S_ij[index_i,index_j] == l && new_S_ij[index_j,index_i] ==m){
              if(Y[index_i, index_j] ==1){
                sum_y_ij <- sum_y_ij+1
              }
            }
          }
        }
        
        #counting number of links among clusters i,j
        sum_y_ij <- 0 
        for(u in 1:length(empty_p)){
          if(Y[empty_p[u],empty_q[u]] ==1){
            sum_y_ij <- sum_y_ij + 1
          }
        }
        
        al1 <- c(1,1,1,1,1)
        bl1 <- c(1,1,1,1,1)
        #updating the C matrix of the block probabilities
        a_l_new <- a + sum_y_ij
        b_l_new <- b + sum_y_ij
        
        a_l <- my_C[which(S_ij[,i,j]>0)] + sum_y_ij
        b_l <- my_C[which(R_ij[,i,j]>0)]+ sum_y_ij
        my_C[which(S_ij[,i,j]>0),which(R_ij[,i,j]>0)] <- rbeta(1,a_l,b_l)
        C_MCMC[,,t] <- my_C
      }
    }
    print(t)
  }
  return(Sij_MCMC)
}