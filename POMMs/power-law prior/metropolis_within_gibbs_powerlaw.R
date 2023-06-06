
alpha_proposal_prime =function(mu, sigma){
  sample_norm_trunc(1,mu,sigma,0.01,10)
}
results_list_powerlaw1 <- list()

alpha_list = c(.5,1,1.5)
K_list= c(3,5,10)

N=100
a=1
b=1
K_max=4
M=4000
max_clust=3
min_clust=3
max_number_games = 10
N_iter = 2000
alpha= 1.5
beta_max = .75


############
set.seed(13)

############
#SIMULATING THE TOURNAMENT
synth_power = simulating_tournament_powerlaw(N=N,alpha = alpha,beta_max = beta_max,min_clust =min_clust ,max_clust = max_clust , M = M, n_ij_max =max_number_games )

#synth = simulating_tournament_test(N=N,alpha = 1,beta = 3,min_clust = min_clust,max_clust = max_clust , M = M, n_ij_max =max_number_games )
synth_matches = synth_power$matches_results
synth_players = synth_power$z_true
synth_p = synth_power$p_true
K_true = synth_power$K_true

test = df_aux_fast(synth_matches,synth_players,synth_p)
n_test = table(synth_players$z)
n_test
cor(test$y_ij/test$n_ij, test$p_ij)
###################

#########
#INITIALIZING THE MCMC
##########
set.seed(23)

labels_available = 1:K_true
#z_proposal = sample(x = labels_available,size = N, replace =  T)
y_ij_matrix= matrix(0,N,N)
for(i in 1:M){
  y_ij_matrix[synth_matches$player_1[i], synth_matches$player_2[i]] =  synth_matches$y_ij[i]
  y_ij_matrix[synth_matches$player_2[i], synth_matches$player_1[i]] =  synth_matches$n_ij[i] - synth_matches$y_ij[i]
}
z_proposal = kmeans(x = y_ij_matrix,centers =K_true)$cluster


z_current = data.frame(id = synth_players$id, z = z_proposal)
df_current = df_aux_fast(synth_matches, z_current, synth_p)
n_current = table(z_current$z)
A_seq_z = matrix(0, 1,N_iter)
z_seq = matrix(0, N, N_iter)
z_seq[, 1] =  z_current$z
A_seq_z[1] = get_A(df_current$n_ij, df_current$y_ij, df_current$p_ij) + dir_multinom_d(N, K_true, n_current)

max_reassign = 3  # maximum number of nodes to reassign
min_acc_ratio = 0.3  # minimum acceptance ratio
reassign_nodes = 2  # initial number of nodes to reassign
acc.count_z=0

##########beta set up



#containers
sigma0=0.1
reassign_nodes_container = matrix(0,1,N_iter*N)
sigma_prime_container = matrix(0,1,N_iter)
alpha_seq= matrix(0,1,N_iter)
A_seq_p = matrix(0,1,N_iter)
p_seq= array(0,c(K_true,K_true,N_iter))
#sampling first value
alpha_current <- proposal_distribution(1,sigma0)
p_current = simulating_POMM_powerlaw(K = K_true,alpha = alpha_current,beta_max = beta_max)
#initializing
data_current = df_aux_fast(synth_matches,synth_players,p_current$matrix)
alpha_seq[1] = alpha_current
A_seq_p[1]= (get_A(data_current$n_ij, data_current$y_ij,data_current$p_ij)+
             l_like_p_ij(p_current$matrix,p_current$truncations) )
sigma_prime_container[1,1]=.1
acc.count_p=0
pb=txtProgressBar(min=1,max=N_iter)
j=2
#######
##METROPOLIS-HASTINGS

while (j < N_iter + 1) {
  
  
  #######z update
  setTxtProgressBar(pb, j)
  z_prime = z_current
  df_prime = df_current
  n_prime = n_current
  for (ii in 1:N) {
    reassign_nodes <- sample(size = 1, x=c(1:max_reassign), p=c(0.5,seq(from = 0.3,to= 0.01,by= ((0.01 - 0.3 )/(max_reassign - 2)))))
    reassign_nodes_container[1,(ii + N*(j-2))] = reassign_nodes
    
    if (reassign_nodes == 1) {
      z_scanning = get_proposal1(z_prime, labels_available)
    } else {
      z_scanning = get_proposal1(z_prime, labels_available)
    }
    df_scanning = df_aux_fast(synth_matches, z_scanning, p_current$matrix)
    n_scanning = table(factor(z_scanning$z, levels = labels_available))
    #r = (get_A(df_scanning$n_ij, df_scanning$y_ij, df_scanning$p_ij) + dir_multinom_d(N, labels_available, n_scanning)) -
      #(get_A(df_prime$n_ij, df_prime$y_ij, df_prime$p_ij) + dir_multinom_d(N, labels_available, n_prime))
    r = (get_A(df_scanning$n_ij, df_scanning$y_ij, df_scanning$p_ij)) -
      (get_A(df_prime$n_ij, df_prime$y_ij, df_prime$p_ij))
    
    u = runif(1)
    if (log(u) < r) {
      acc.count_z = acc.count_z + 1
      z_prime = z_scanning
      df_prime = df_scanning
      n_prime = n_scanning
    }
    # } else {
    #   
    #   # if (acc.count / j < min_acc_ratio && reassign_nodes < max_reassign) {
    #   #   reassign_nodes = reassign_nodes + 1
    #   # }
    # }
    A_seq_z[1,j] = get_A(df_prime$n_ij, df_prime$y_ij, df_prime$p_ij) + dir_multinom_d(N, labels_available, n_prime)
    z_current = z_prime
    df_current = df_prime
    n_current = n_prime
    z_seq[, j-1] = z_prime$z
  }
  
#######beta update
  #updating sigma
  #sigma_prime= sample_norm_trunc(1, sigma_prime_container[1,j-1], s = (acc.count_p+1)/j,a = 0.001,b=1)
  sigma_prime = runif(1,0.01,0.3)
  sigma_prime_container[1,j] = sigma_prime
  #proposing a new alpha
  alpha_prime <- alpha_proposal_prime(alpha_current,sigma =sigma_prime)
  #generating a proposal matrix
  p_prime = simulating_POMM_powerlaw(K_true,alpha_prime,beta_max)
  #mapping it into the aux dataframe
  data_prime = df_aux_fast(synth_matches,synth_players,p_prime$matrix)
  
  #evaluating the likelihood
  r = (get_A(data_prime$n_ij, data_prime$y_ij,data_prime$p_ij)+
         l_like_p_ij(p_prime$matrix,p_prime$truncations)) -
    (get_A(data_current$n_ij, data_current$y_ij,data_current$p_ij)+
       l_like_p_ij(p_current$matrix,p_current$truncations))

  u=runif(1)
  if(log(u)<r){
    acc.count_p= acc.count_p+1
    data_current = data_prime
    p_current = p_prime
    alpha_current = alpha_prime
    alpha_seq[j] = alpha_prime
    p_seq[,,j] = p_prime$matrix
    A_seq_p[j] =  (get_A(data_prime$n_ij, data_prime$y_ij,data_prime$p_ij)+
                   l_like_p_ij(p_prime$matrix,p_prime$truncations))
  }else{
    alpha_seq[j] = alpha_current
    p_seq[,,j] = p_current$matrix
    A_seq_p[1,j] = (get_A(data_current$n_ij, data_current$y_ij,data_current$p_ij)+
                    l_like_p_ij(p_current$matrix,p_current$truncations))
  }
  
  j=j+1
}

iteration_title = paste("K",K_true, "alpha" ,alpha)
iteration_title= gsub(" ", "",as.character(iteration_title))
iteration_title= gsub(".", "_",as.character(iteration_title),fixed = T)

similarity<- pr_cc(z_seq[,-c(1:(N_iter*0.5))])


#POINT ESTIMATES


#POINT ESTIMATES
MAP = which(A_seq_z[1:(N_iter*0.7-1)] == max(A_seq_z[1:(N_iter*0.7-1)]))


adj.rand = adj.rand.index(synth_players$z,z_seq[,MAP[1]])
diag(similarity)=1

adj.rand.VI= adj.rand.index(synth_players$z,minVI(similarity,method = "avg")$cl)

similarity_plot(similarity,z_0 = synth_players$z, z_est =synth_players$z)


acc.percentage_z = acc.count_z/(N_iter*N) *100

acc.percentage_p = acc.count_p/(N_iter) *100

ts.plot(A_seq_z[-c(1:(N_iter*0.5))])

acf(A_seq_z[-c(1:(N_iter*0.9))])


#checking p

#checking the mixing of the distribution
ts.plot(A_seq_p[-c(1:(N_iter*0.5))])

plot(ts(alpha_seq[-c(1:(N_iter*0.2))]))
abline(h = alpha, col = "red", lty = 2)

#looking at the MSE between the MAP and the true one
mse_table = matrix(0,K_true,K_true)
for(i in 1:K_true){
  for(j in 1:K_true){
    mse_table[i,j] = (mean(p_seq[i,j,]) - synth_p[i,j])**2
  }}
library(pander)
mse_table%>%pander()

# 
# burnin_p = p_seq[,,20001:30000]
results_list_powerlaw1[[iteration_title]] <- list(
  alpha_seq = alpha_seq,
  acc.count_p = acc.count_p,
  acc.count_z = acc.count_z,
  A_seq_p = A_seq_p,
  A_seq_z = A_seq_z,
  p_seq_power = p_seq,
  z_seq_power = z_seq,
  true_z = synth_players,
  similarity_matrix = similarity,
  adj_rand_index = adj.rand,
  mse_table = mse_table)



sqrt(sum((p_seq[,,which(A_seq_p==max(A_seq_p))[1]] - synth_p)**2))

mean(p_seq[1,3,-(N_iter*0.5)])

burnin_p = p_seq[,,-(N_iter*0.5)]

plots = list()
for(i in 1:K_true) {
  for(j in 1:K_true) {
    y_try = data.frame(y = as.vector(burnin_p[i, j,]))
    p1 = ggplot(y_try, aes(y)) +
      geom_density(fill = "dodgerblue", alpha = 0.5) +
      scale_x_log10() +
      geom_vline(xintercept = synth_p[i, j], color = "red")
    plots[[length(plots) + 1]] <- p1
  }
}
p_combined = wrap_plots(plots, ncol = 3, nrow = 3)
p_combined

plots[[2]]
results_list_powerlaw1$K10alpha0_1$adj_rand_index

