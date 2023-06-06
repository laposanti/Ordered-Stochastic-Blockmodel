


mse_iterations_pomm = list()
alpha_iterations_pomm = list()
rand_index_list_pomm = list()



setwd("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/comparison with simplified model/new_study")

source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")

library(dplyr)
library(fossil)
library(ggplot2)
library(collpcm)

############
#SIMULATING THE TOURNAMENT
############
alpha_list = c(.5,1,1.5)
K_list = c(3,5,8)
my_seed =13
while(my_seed<=23){
# for(p in 1:length(alpha_list)){
#   for(pp in 1:length(K_list)){

    #Hyperparameters
    N=100
    a=1
    b=1
    K_max=4
    M=3800
    max_clust=3
    min_clust= 3
    max_number_games = 100
    N_iter = 10000
    alpha= 0.5
    gamma_vec = rep(1,min_clust)
    beta_max = .8
    
    set.seed(my_seed) 
    
    iteration_title = paste("seed", my_seed, "K",min_clust, "alpha" ,alpha, sep = "")
    
    print(iteration_title)


synth_power = simulating_tournament_powerlaw(N=N,alpha = alpha,beta_max = beta_max,min_clust = min_clust,max_clust = max_clust , M = M, n_ij_max =max_number_games, gamma_vec = gamma_vec )
#synth = simulating_tournament_test(N=N,alpha = 1,beta = 3,min_clust = min_clust,max_clust = max_clust , M = M, n_ij_max =max_number_games )
synth_matches = synth_power$matches_results
synth_players = synth_power$z_true
synth_p = synth_power$p_true
K_true = synth_power$K_true


df_fast =df_aux_fast(synth_matches,synth_players,synth_p)


p_ij.true = synth_p
z.true = synth_players$z
z_mat_true=vec2mat(z.true)
matrix_z_p_true = p_ij.true%*%t(z_mat_true)
p_n_true = z_mat_true%*%matrix_z_p_true


n_ij_matrix = matrix(0,N,N)
for(i in 1:M){
  n_ij_matrix[synth_matches$player_1[i], synth_matches$player_2[i]] =  synth_matches$n_ij[i]
  n_ij_matrix[synth_matches$player_2[i], synth_matches$player_1[i]] =  synth_matches$n_ij[i]
}

y_ij_matrix= matrix(0,N,N)
for(i in 1:M){
  y_ij_matrix[synth_matches$player_1[i], synth_matches$player_2[i]] =  synth_matches$y_ij[i]
  y_ij_matrix[synth_matches$player_2[i], synth_matches$player_1[i]] =  synth_matches$n_ij[i] - synth_matches$y_ij[i]
}




z.true = synth_players$z
labels_available=1:K_true




#selecting just those values such that there is inside a non-0 entry
upper.tri_n_ij = upper.tri(n_ij_matrix)
upper.tri_y_ij = upper.tri(y_ij_matrix)
non_negative_n_ij = which(upper.tri_n_ij & n_ij_matrix > 0, arr.ind = T)
non_negative_y_ij = which(upper.tri_y_ij & y_ij_matrix > 0, arr.ind = T)


#retrieving the non-0 values
n_ij = n_ij_matrix[non_negative_n_ij]
y_ij = y_ij_matrix[non_negative_n_ij]


p_ij.true = synth_p
z.true = synth_players$z
z_mat_true=vec2mat(z.true)
matrix_z_p_true = p_ij.true%*%t(z_mat_true)
p_n_true = z_mat_true%*%matrix_z_p_true
p_ij_true =  p_n_true[non_negative_n_ij]

if(round(sum(dbinom(y_ij,n_ij,p_ij_true,log=T)),2) == round(get_A(df_fast$n_ij,df_fast$y_ij,pij = df_fast$p_ij),2)){
  print("METROPOLIS: CLEAR!")}



#Containers for z
#---
z.container = matrix(0,N,N_iter)
A_seq = matrix(0,N_iter,1)

#Initialization for the z vector
#---
init = kmeans(x = y_ij_matrix,centers =K_true)$cluster
adj.rand.index(init,z.true)
z_current= init
z.container[,1] = z_current

#Containers for p
#---
alpha.container = matrix(0,N_iter,1)
p.container = array(0,c(K_true,K_true,N_iter))


#Initialization for the p matrix
#---
alpha_current = 1
p_current = simulating_POMM_powerlaw(K_true,alpha_current,beta_max = beta_max)

alpha.container[1] = alpha_current
p.container[,,1] = p_current$matrix



#Setting up quantities needed within computations

# here we are transforming a Nx1 vector, containg labels 1...K into
# NXK matrix. z_mat_current[i,k] =1 if node i is in cluster k
z_mat_current= vec2mat(z_current)
n_k_current = colSums(z_mat_current)
#creating an NxN matrix where p_n_current[i,j] is the probability that player i wins vs player j
matrix_z_p_current = p_current$matrix%*%t(z_mat_current)
p_n_current = z_mat_current%*%matrix_z_p_current
p_ij_current =  p_n_current[non_negative_n_ij]


A_seq[1] = sum(dbinom(y_ij, n_ij, p_ij_current, log=T)) + l_like_p_ij(p_current$matrix,p_current$truncations)



#containers for the counts of accepted proposals
acc.count_z = 0
acc.count_p = 0

#setting time tracker
pb=txtProgressBar(min=1,max=N_iter)
j=2

#READY TO BOMB!
while (j < N_iter + 1) {
  setTxtProgressBar(pb, j)
  
  #Complete sweeep of z vector
  #----
  sweeping_order = sample(x=c(1:N),size=N, replace=F)
  for(i in sweeping_order){
    #proposing a new label
    new_label = sample(x=labels_available,size=1)
    
    #updating labels
    z_prime = z_current
    z_prime[i] = new_label
    
    #computing the new victory probabilities
    while(TRUE){
      z_mat_prime= vec2mat(z_prime)
      if(ncol(z_mat_prime) == K_true){
        break
      } else {
        #if there is an error, resample new_label and try again
        new_label = sample(x=labels_available,size=1)
        z_prime[i] = new_label
      }
    }
    matrix_z_prime_p_current = p_current$matrix%*%t(z_mat_prime)
    p_n_z_prime = z_mat_prime%*%matrix_z_prime_p_current
    p_ij_z_prime =  p_n_z_prime[non_negative_n_ij]
    
    n_k_prime = colSums(z_mat_prime)
    #acceptance ratio
    r = (sum(dbinom(y_ij, n_ij, p_ij_z_prime, log=T)) + ddirichlet_multinomial(N,K_true,n_k = n_k_prime,my_alpha = gamma_vec)) - 
      (sum(dbinom(y_ij, n_ij, p_ij_current, log = T)) + ddirichlet_multinomial(N,K_true,n_k = n_k_current ,my_alpha = gamma_vec))
    
    alpha_r = min(1, exp(r))
    u = runif(1)
    
    #if accepted
    if(u<alpha_r){
      acc.count_z = acc.count_z + 1
      z_current=z_prime
      n_k_current = n_k_prime
      p_ij_current= p_ij_z_prime
      #if not accepted
    }
   
  }
  z.container[, j] = z_current
  #Update of the P matrix
  #----
  sigma_prime = .2
  
  #proposing a new alpha
  alpha_prime <- sample_norm_trunc(1,m = alpha_current,s =sigma_prime,a = 0.005,b=3)

  #generating a proposal matrix
  p_prime = simulating_POMM_powerlaw(K_true,alpha_prime,beta_max)
  
  

  #Updating p_ij_prime with the last membership
  z_mat_current = vec2mat(z_current)
  matrix_z_p_prime  = p_prime$matrix%*%t(z_mat_current)
  p_n_prime = z_mat_current%*%matrix_z_p_prime
  p_ij_prime =  p_n_prime[non_negative_n_ij]
  
  r = (sum(dbinom(y_ij, n_ij, p_ij_prime, log=T)) +  l_like_p_ij(p_prime$matrix,p_prime$truncations)) - 
    (sum(dbinom(y_ij, n_ij, p_ij_current, log = T)) + l_like_p_ij(p_current$matrix,p_current$truncations) )
  
  alpha_r = min(1, exp(r))
  u = runif(1)
  if(u<alpha_r){
    #counting number of accepted proposals
    acc.count_p = acc.count_p+1
    #updating quantities
    p_ij_current = p_ij_prime
    alpha_current = alpha_prime
    p_current = p_prime
  }
  #storing results for diagnostics
  p.container[,,j] = p_current$matrix
  alpha.container[j]= alpha_current
  
  A_seq[j] = sum(dbinom(y_ij, n_ij, p_ij_current, log=T)) +  l_like_p_ij(p_current$matrix,p_current$truncations)
  j=j+1
}


#diagnostic for z
#-----
#proposal performance
acceptance_rate= acc.count_z/(N*N_iter)*100

#mixing
# ts.plot(A_seq[-c(1:N_iter*0.5)], main="Traceplot of p(y_ij|z,P) p(z) p(P)",
#         xlab = "Iterations", ylab ="p(y_ij|z,P) p(z) p(P)")
# dev.off()
# acf(A_seq[-c(1:N_iter*0.5)],main="Autocorrelation plot",
#     xlab = "Lag", ylab ="ACF")
# dev.off()


#estimates
similarity_matrix = pr_cc(z.container[,-c(1:N_iter*0.5)])

point_est = minVI(similarity_matrix)$cl
adj.rand.index(point_est, z.true)


similarity_plot(y_ij_matrix,z.true,z.true)
dev.off()


similarity_plot(similarity_matrix,z.true,z.true)
dev.off()

#diagnostic for P
#-----

acceptance_rate_p  = (acc.count_p/N_iter)*100

plot(ts(alpha.container), main = "Traceplot of alpha values")
abline(h = .5, col = "red", lty = 2)



# Create a data frame containing the observations
df_alpha <- data.frame(alpha_est = alpha.container)

# Create the plot
# ggplot(df_alpha, aes(x = alpha_est)) +
#   geom_density(fill = "purple", alpha = 0.5) +
#   geom_vline(aes(xintercept = mean(alpha_est)), color = "blue") +
#   geom_vline(aes(xintercept = alpha), color = "red")

for(i in 1:K_true){
  for(j in 1:K_true){
    print(mean(p.container[i,j,]) - synth_p[i,j])
  }
}


mse_table = matrix(0,K_true,K_true)
for(i in 1:K_true){
  for(j in 1:K_true){
    mse_table[i,j]= (mean(p.container[i,j,]) - synth_p[i,j])
  }
}

mse_table%>% pander::pander()

mse_total = sum(mse_table[upper.tri(mse_table)]**2)/((K_true*(K_true-1))/2)


        
burnin_p = p.container[,,-(N_iter*0.5)]

# plots = list()
# for(i in 1:K_true) {
#   for(j in 1:K_true) {
#     y_try = data.frame(y = as.vector(burnin_p[i, j,]))
#     p1 = ggplot(y_try, aes(y)) +
#       geom_density(fill = "dodgerblue", alpha = 0.5) +
#       scale_x_log10() +
#       geom_vline(xintercept = synth_p[i, j], color = "red")+
#       xlab("probability") +
#       ylab("Density") +
#       ggtitle(paste("Density plot of entry ", i, ",", j, sep = ""))
#     
#     plots[[length(plots) + 1]] <- p1
#   }
# }
# p_combined = patchwork::wrap_plots(plots, ncol = 5, nrow = 5)
# p_combined


png(as.character(paste("adjacencypomm",iteration_title,".png",sep = "")),width = 800, height = 540)
similarity_plot(y_ij_matrix,z.true,z.true)
dev.off()

 png(as.character(paste("similaritypomm",iteration_title,".png",sep = "")),width = 800, height = 540)
similarity_plot(similarity_matrix,z.true,z.true)
 dev.off()

 # png(as.character(paste("heatmap_mse_pomm",iteration_title,".png",sep = "")),width = 400, height = 400)
 # heat_map_blue(matrix = mse_table, "heatmap_mse")
 # dev.off()


 mse_iterations_pomm[iteration_title] = sqrt(mse_total)
 rand_index_list_pomm[iteration_title] = adj.rand.index(point_est, z.true)
 alpha_iterations_pomm[iteration_title] = mean(alpha.container[-c(1:N_iter*0.5)])
 my_seed = my_seed+1
 
 
  }




boxplot(as.vector(unlist(rand_index_list_pomm)),
        main= "Adj.Rand Index for 10 samples, K=8, alpha=0.5",
        ylab = "Adj.Rand.Index")

mse_vec = as.vector(unlist(mse_iterations_pomm))
boxplot(mse_vec,main= "RMSE for 10 samples, K=8, alpha=0.5",
        ylab = "RMSE")
alpha_error = as.vector(unlist(alpha_iterations_pomm)) - alpha

boxplot(alpha_error,main= "Alpha error for 10 samples, K=8, alpha=0.5",
        ylab = "alpha hat -alpha")



k8_adj_rand_index = unlist(rand_index_list_pomm[c(1:11)]) 
k5_adj_rand_index = unlist(rand_index_list_pomm[c(12:22)])
k3_adj_rand_index = unlist(rand_index_list_pomm[c(23:33)])

k8_rmse_P = unlist(mse_iterations_pomm[c(1:11)])
k5_rmse_P = unlist(mse_iterations_pomm[c(12:22)])
k3_rmse_P = unlist(mse_iterations_pomm[c(23:33)])

k8_rmse_alpha = unlist(alpha_iterations_pomm[c(1:11)])- alpha
k5_rmse_alpha = unlist(alpha_iterations_pomm[c(12:22)])- alpha
k3_rmse_alpha = unlist(alpha_iterations_pomm[c(23:33)])- alpha
mean(k3_rmse_alpha)

mean(k8_rmse_alpha)
library(ggplot2)
library(dplyr)
library(tidyr)

#####
#PLOT RESULTS OF THE SIMULATIONS
######


# Combine the data into a data frame
data_rand <- data.frame(state = c(rep("K=3", 11), rep("K=5", 11), rep("K=8", 11)),
                   indicator = c(k3_adj_rand_index, k5_adj_rand_index, k8_adj_rand_index))

# Calculate the mean and standard deviation of the indicator for each state
means_rand <- data_rand %>%
  group_by(state) %>%
  summarize(mean = median(indicator),
            sd = sd(indicator),
            se = sd / sqrt(n()))

# Plot the means and confidence intervals
ggplot(means_rand, aes(x = state, y = mean))  +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean *0 +1),
                width = 0.4, color = c("blue","red","purple")) +
  geom_point(aes(y = mean), color = "red")+
  labs(title = "Mean and 95% CI of Adjusted Rand Index by K",
       x = "State",
       y = "Mean Indicator") +
  theme_bw()

# Combine the data into a data frame
data_RMSE_P <- data.frame(state = c(rep("K=3", 11), rep("K=5", 11), rep("K=8", 11)),
                        indicator = c(k3_rmse_P, k5_rmse_P, k8_rmse_P))

# Calculate the mean and standard deviation of the indicator for each state
means_RMSE_P <- data_RMSE_P %>%
  group_by(state) %>%
  summarize(mean = median(indicator),
            sd = sd(indicator),
            se = sd / sqrt(n()))

# Plot the means and confidence intervals
ggplot(means_RMSE_P, aes(x = state, y = mean))  +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se),
                width = 0.4, color = c("blue","red","purple")) +
  geom_point(aes(y = mean), color = "red")+
  labs(title = "Mean and 95% CI of RMSE for P by K",
       x = "State",
       y = "Mean Indicator") +
  theme_bw()


# Combine the data into a data frame
data_RMSE_alpha <- data.frame(state = c(rep("K=3", 11), rep("K=5", 11), rep("K=8", 11)),
                          indicator = c(k3_rmse_alpha, k5_rmse_alpha, k8_rmse_alpha))

# Calculate the mean and standard deviation of the indicator for each state
data_RMSE_alpha <- data_RMSE_alpha %>%
  group_by(state) %>%
  summarize(mean = median(indicator),
            sd = sd(indicator),
            se = sd / sqrt(n()))

# Plot the means and confidence intervals
ggplot(data_RMSE_alpha, aes(x = state, y = mean))  +
  geom_errorbar(aes(ymin = mean*0, ymax = mean + 1.96 * se),
                width = 0.4, color = c("blue","red","purple")) +
  geom_point(aes(y = mean), color = "red")+
  labs(title = "Mean and 95% CI of RMSE_alpha by K",
       x = "State",
       y = "Mean Indicator") +
  theme_bw()
