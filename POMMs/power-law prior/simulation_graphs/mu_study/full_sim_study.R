
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(fossil)
############
#simulation study
############


results_list_powerlaw <- list()
results_list_simple <- list()
alpha_list = c(.5,1,1.5)
K_list= c(3,5,10)
for(t in 1:length(alpha_list)){
  for(tt in 1:length(K_list)){
    
    ############
    #SIMULATING THE TOURNAMENT
    ############
    
    N=100
    a=1
    b=1
    K_max=4
    M=3900
    max_clust= K_list[tt]
    min_clust= K_list[tt]
    max_number_games = 50
    N_iter = 10000
    alpha= alpha_list[t]
    beta_max = .75
    ############
    set.seed(13)
    ############
    #SIMULATING THE TOURNAMENT
    synth_power = simulating_tournament_powerlaw(N=N,alpha = alpha,beta_max = beta_max,min_clust = min_clust,max_clust = max_clust , M = M, n_ij_max =max_number_games )
    #synth = simulating_tournament_test(N=N,alpha = 1,beta = 3,min_clust = min_clust,max_clust = max_clust , M = M, n_ij_max =max_number_games )
    synth_matches = synth_power$matches_results
    synth_players = synth_power$z_true
    synth_p = synth_power$p_true
    K_true = synth_power$K_true
    
    test = df_aux_fast(synth_matches,synth_players,synth_p)
    n_test = table(synth_players$z)
    n_test
    cor(test$y_ij/test$n_ij, test$p_ij)
    
    
    #############
    #simple model
    #########
    
    
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
    
    
    iteration_title = paste("K",K_true, "alpha" ,alpha)
    iteration_title= gsub(" ", "",as.character(iteration_title))
    iteration_title= gsub(".", "_",as.character(iteration_title),fixed = T)
    iteration_title
    
    similarity_matrix<- pr_cc(z_seq[,-c(1:(N_iter*0.25))])
    
    
    #POINT ESTIMATES
    MAP = which(A_seq[1:(N_iter-1)] == max(A_seq[1:(N_iter-1)]))
    
    adj_rand_index = adj.rand.index(synth_players$z,z_seq[,MAP[1]])
    
    similarity_plot(similarity_matrix,z_0 = synth_players$z, z_est =synth_players$z)
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
    print("simple")
    mse_table%>%pander()
    
    
    
    results_list_simple[[iteration_title]] <- list(
      acc.count_p = acc.count_p,
      acc.count_z = acc.count_z,
      A_seq_p = A_seq_p,
      A_seq_z = A_seq_z,
      p_seq = p_seq,
      z_seq = z_seq,
      true_z = synth_players,
      similarity_matrix = similarity_matrix,
      adj_rand_index = adj_rand_index,
      mse_table = mse_table)
    
    
    
    
    
    
    
    ###############
    #power law model
    ###########
    
  
    
    
    #########
    #INITIALIZING THE MCMC
    ##########
    
    
    labels_available = 1:K_true
    z_proposal = sample(x = labels_available,size = N, replace =  T)
    
    
    z_current = data.frame(id = synth_players$id, z = z_proposal)
    df_current = df_aux_fast(synth_matches, z_current, synth_p)
    n_current = table(z_current$z)
    A_seq_z = matrix(0, N_iter, 1)
    z_seq_power = matrix(0, N, N_iter)
    z_seq_power[, 1] =  z_current$z
    A_seq_z[1] = get_A(df_current$n_ij, df_current$y_ij, df_current$p_ij) + dir_multinom_d(N, K_true, n_current)
    
    max_reassign = 3  # maximum number of nodes to reassign
    min_acc_ratio = 0.3  # minimum acceptance ratio
    reassign_nodes = 2  # initial number of nodes to reassign
    acc.count_z=0
    
    ##########beta set up
    
    
    
    #containers
    alpha_seq= matrix(0,1,N_iter)
    A_seq_p = matrix(0,1,N_iter)
    p_seq_power= array(0,c(K_true,K_true,N_iter))
    #sampling first value
    alpha_current <- proposal_distribution(1,sigma0)
    p_current = simulating_POMM_powerlaw(K = K_true,alpha = alpha_current,beta_max = beta_max)
    #initializing
    data_current = df_aux_fast(synth_matches,synth_players,p_current$matrix)
    alpha_seq[1] = alpha_current
    A_seq_p[1]= (get_A(data_current$n_ij, data_current$y_ij,data_current$p_ij)+
                 l_like_p_ij(p_current$matrix,p_current$truncations) )
    
    acc.count_p=0
    pb=txtProgressBar(min=1,max=N_iter)
    j=2
    #######
    ##METROPOLIS-HASTINGS
    
    while (j < N_iter + 1) {
      setTxtProgressBar(pb, j)
      #######z update
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
        r = (get_A(df_scanning$n_ij, df_scanning$y_ij, df_scanning$p_ij) + dir_multinom_d(N, labels_available, n_scanning)) -
          (get_A(df_prime$n_ij, df_prime$y_ij, df_prime$p_ij) + dir_multinom_d(N, labels_available, n_prime))
        u = runif(1)
        if (log(u) < r) {
          acc.count_z = acc.count_z + 1
          z_prime = z_scanning
          df_prime = df_scanning
          n_prime = n_scanning
          A_seq_z[j] = get_A(df_prime$n_ij, df_prime$y_ij, df_prime$p_ij) + dir_multinom_d(N, labels_available, n_prime)
        } else {
          A_seq_z[j] = A_seq_z[j - 1]
          # if (acc.count / j < min_acc_ratio && reassign_nodes < max_reassign) {
          #   reassign_nodes = reassign_nodes + 1
          # }
        }
        z_current = z_prime
        df_current = df_prime
        n_current = n_prime
        z_seq_power[, j - 1] = z_prime$z
      }
      
      #######beta update
      #updating sigma
      sigma_prime= runif(1,0,2)
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
        p_seq_power[,,j] = p_prime$matrix
        A_seq_p[j] =  (get_A(data_prime$n_ij, data_prime$y_ij,data_prime$p_ij)+
                         l_like_p_ij(p_prime$matrix,p_prime$truncations))
      }else{
        alpha_seq[j] = alpha_seq[j-1]
        p_seq_power[,,j] = p_current$matrix
        A_seq_p[j] = A_seq_p[j-1]
      }
      
      j=j+1
    }
    
    iteration_title = paste("K",K_true, "alpha" ,alpha)
    iteration_title= gsub(" ", "",as.character(iteration_title))
    iteration_title= gsub(".", "_",as.character(iteration_title),fixed = T)
    #checking z
    
    similarity_matrix<- pr_cc(z_seq_power[,-c(1:(N_iter*0.25))])
    

    #POINT ESTIMATES
    
    
    MAP = which(A_seq_z[1:(N_iter-1)] == max(A_seq_z[1:(N_iter-1)]))
    
    adj_rand_index = adj.rand.index(synth_players$z,z_seq_power[,MAP[1]])
    
    
    similarity_plot(similarity_matrix,z_0 = synth_players$z, z_est =synth_players$z)
    ts.plot(A_seq_z[-c(1:(N_iter*0.25))])
    #checking p
    
    #checking the mixing of the distribution
    ts.plot(A_seq_p[-c(1:(N_iter*0.25))])
    
    #looking at the MSE between the MAP and the true one
    mse_table = matrix(0,K_true,K_true)
    for(i in 1:K_true){
      for(j in 1:K_true){
        mse_table[i,j] = (mean(p_seq_power[i,j,]) - synth_p[i,j])**2
      }}
    library(pander)
    mse_table%>%pander()
    
  
    
    results_list_powerlaw[[iteration_title]] <- list(
      alpha_seq = alpha_seq,
      acc.count_p = acc.count_p,
      acc.count_z = acc.count_z,
      A_seq_p = A_seq_p,
      A_seq_z = A_seq_z,
      p_seq_power = p_seq,
      z_seq_power = z_seq,
      true_z = synth_players,
      similarity_matrix = similarity_matrix,
      adj_rand_index = adj_rand_index,
      mse_table = mse_table)
  }}



# Define a function that takes a matrix and a membership vector and returns a heatmap
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

diag(results_list_powerlaw$K5alpha1_5$similarity_matrix)=1
isSymmetric(results_list_powerlaw$K5alpha1_5$similarity_matrix)

z_est= minVI(results_list_powerlaw$K5alpha1_5$similarity_matrix, method = "avg", max.k = 20)$cl

# Example usage with random data
set.seed(123)
mat <- matrix(rnorm(100), nrow = 10)
memb <- sample(1:3, 10, replace = TRUE)

similarity_plot(results_list_simple$K3alpha1$similarity_matrix,results_list_simple$K3alpha1$true_z$z,results_list_simple$K3alpha1$true_z$z)


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




sim1 = results_list_powerlaw$K3alpha1$similarity_matrix
true1 = results_list_powerlaw$K3alpha1$true_z$z
rand1 = results_list_powerlaw$K3alpha1$adj_rand_index

sim2 = results_list_simple$K5alpha1$similarity_matrix
true2 = results_list_simple$K5alpha1$true_z$z
rand2 = results_list_simple$K3alpha1$adj_rand_index

similarity_plot(sim1,true1,true1)
similarity_plot(sim2,true2,true2)





