

source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")

p_proposal = function(p, sigma_p,K){
  p_new <- matrix(sample_norm_trunc(K*K, p, sigma_p,0,1),K,K)
  return(p_new)
}
############
#SIMULATING THE TOURNAMENT
############

#Hyperparameters
N=100
a=1
b=1
K_max=4
M=3389
max_clust=3
min_clust=3
max_number_games = 4
N_iter = 10000
alpha=1
beta_max = .8
gamma_vec = rep(1,min_clust)

set.seed(123)
synth_power = simulating_tournament_powerlaw(N=N,alpha = alpha,beta_max = beta_max,min_clust = min_clust,max_clust = max_clust , M = M, n_ij_max =max_number_games, gamma_vec = gamma_vec )
#synth = simulating_tournament_test(N=N,alpha = 1,beta = 3,min_clust = min_clust,max_clust = max_clust , M = M, n_ij_max =max_number_games )
synth_matches = synth_power$matches_results
synth_players = synth_power$z_true
synth_p = synth_power$p_true
K_true = synth_power$K_true




###########
#Containers
###########

p_seq= array(0,c(K_true,K_true,N_iter))
z_seq = matrix(0, N, N_iter)

A_seq_p = matrix(0,N_iter,1)
A_seq_z = matrix(0, N_iter, 1)
#############
#Initializing quantities
###########


p_current <- p_proposal(rbeta(K_true**2,1,1),2,K = K_true)
p_seq[,,1] = p_current

labels_available = 1:K_true
z_proposal = sample(x = labels_available,size = N, replace =  T)
z_current = data.frame(id = synth_players$id, z = z_proposal)

z_seq[, 1] =  z_current$z

A_seq_p[1]= get_A(df_current$n_ij, df_current$y_ij,df_current$p_ij)+get_B(p_current, beta_0)
A_seq_z[1] = get_A(df_current$n_ij, df_current$y_ij, df_current$p_ij) + dir_multinom_d(N, K_true, n_current)


df_current = df_aux_fast(synth_matches, z_current, p_current)
n_current = table(z_current$z)




acc.count_z = 0
acc.count_p = 0

j = 2
max_reassign = 3  # maximum number of nodes to reassign
min_acc_ratio = 0.3  # minimum acceptance ratio
reassign_nodes = 2  # initial number of nodes to reassign
pb = txtProgressBar(min = 1, max = N_iter)

#creating an updatable version of z_current, df_current, n_current
z_prime = z_current
df_prime = df_current
n_prime = n_current

while(j<N_iter+1){
  setTxtProgressBar(pb, j)
  ###############
  # z update
  ###############
    
    #we are selecting the number of nodes to reassign
    for (ii in 1:N) {
      reassign_nodes <- sample(size = 1, x=c(1:max_reassign), p=c(0.5,seq(from = 0.3,to= 0.01,by= ((0.01 - 0.3 )/(max_reassign - 2)))))
      if (reassign_nodes == 1) {
        #if just one node is to be reassigned use: get_proposal1
        z_scanning = get_proposal1(z_prime, labels_available)
      } else {
        #if just two or more nodes are to be reassigned use: get_proposal2
        z_scanning = get_proposal2(z_prime, labels_available, reassign_nodes)
      }
      #a container with the newly sampled assignment
      df_scanning = df_aux_fast(synth_matches, z_scanning, p_current)
      n_scanning = table(factor(z_scanning$z, levels = labels_available)) #update also the number of nodes for each block
      #evaluating posterior ratios
      r = (get_A(df_scanning$n_ij, df_scanning$y_ij, df_scanning$p_ij) + dir_multinom_d(N, labels_available, n_scanning)) -
        (get_A(df_prime$n_ij, df_prime$y_ij, df_prime$p_ij) + dir_multinom_d(N, labels_available, n_prime))
      u = runif(1)
      #if proposal gets accepted
      if (log(u) < r) {
        acc.count_z = acc.count_z + 1 # accepted counts + 1
        z_prime = z_scanning # update z_prime
        df_prime = df_scanning #update df_prime
        n_prime = n_scanning #update n_prime
         #storing likelihood value
      } else {
        #does not change likelihood
        # if (acc.count / j < min_acc_ratio && reassign_nodes < max_reassign) {
        #   reassign_nodes = reassign_nodes + 1
        # }
      }
      # at the end of the sweep, store in z_current the fully updated z_prime, df_prime, n_prime
      z_current = z_prime 
      df_current = df_prime
      n_current = n_prime
      z_seq[, j - 1] = z_prime$z
      A_seq_z[j] = get_A(df_current$n_ij, df_current$y_ij, df_current$p_ij) + dir_multinom_d(N, labels_available, n_current)
    }
  
  
  
  
  ######
  #Beta update
  ####

  #proposing a new p
  p_prime = p_proposal(p = p_current,sigma_p = sigma0*exp((acc.count_p+1)/j),K = K_true)
  p_prime[is.na(p_prime)] <- 0
  p_prime = (1 -t(p_prime*upper.tri(p_prime))) * (lower.tri(p_prime, diag = F)*1) + upper.tri(p_prime, diag = T) *p_prime
  diag(p_prime)=0.5
  #mapping it into the aux dataframe
  df_prime = df_aux_fast(synth_matches,z_current,p_prime)
  
  #evaluating the likelihood
  r = (get_A(df_prime$n_ij, df_prime$y_ij,df_prime$p_ij)+get_B(p_prime, beta_0))-
    (get_A(df_current$n_ij, df_current$y_ij,df_current$p_ij)+get_B(p_current, beta_0))
  
  u=runif(1)
  if(log(u)<r){
    acc.count_p=acc.count_p+1
    df_current = df_prime
    p_current = p_prime
    p_seq[,,j] = p_prime
    A_seq_p[j] = get_A(df_prime$n_ij, df_prime$y_ij,df_prime$p_ij)+get_B(p_prime, beta_0)
  }else{
    p_seq[,,j] = p_current
    A_seq_p[j] = get_A(df_current$n_ij, df_current$y_ij,df_current$p_ij)+get_B(p_current, beta_0)
  }
  j=j+1
}

#checking z

similarity<- pr_cc(z_seq[,-c(1:(N_iter*0.25))])


#POINT ESTIMATES


adj.rand.index(synth_players$z,z_seq[,MAP[1]])

similarity_plot(similarity,z_0 = synth_players$z, z_est =synth_players$z)
ts.plot(A_seq_z[-c(1:(N_iter*0.25))])
#checking p

#checking the mixing of the distribution
ts.plot(A_seq_p[-c(1:(N_iter*0.25))])

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




sqrt(sum((p_seq[,,which(A_seq_p==max(A_seq_p))[1]] - synth_p)**2))


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



