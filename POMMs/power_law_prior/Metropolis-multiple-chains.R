

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
#INITIALIZING THE MCMC 4 CHAINS
##########

N.chains= 4
labels.available = 1:K_true
z.samples = array(0, N, c(N.iter,N.chains))
A_seq_z = array(0, 1,c(N.iter,N.chains))

for(i in 1:N.chains){
  z.samples[,1,i] = sample(x = labels_available,size = N, replace =  T)
}

z.current = list()
for(i in 1:N.chains){
  z.current[[i]]= data.frame(id = synth_players$id, z = z.samples[,1,i])
}

df.current = list()
for(i in 1:N.chains){
  df_current[[i]]= df_aux_fast(synth_matches,  z.current[[i]], synth_p)
}
  
n.current = list() 
for(i in 1:N.chains){
  n.current[[i]]= table(z.current[[i]]$z)
}

for(i in 1:N.chains){
  A_seq_z[1,i] = get_A(df_current[[i]]$n_ij, df_current[[i]]$y_ij, df_current[[i]]$p_ij) + dir_multinom_d(N, K_true, n_current[[i]])
}


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