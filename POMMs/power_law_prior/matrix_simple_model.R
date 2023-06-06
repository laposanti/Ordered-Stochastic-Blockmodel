
setwd("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/comparison with simplified model")

mse_iterations = list()
alpha_iterations = list()
rand_index_list = list()

source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
library(dplyr)
library(fossil)
library(ggplot2)
############
#SIMULATING THE TOURNAMENT
############

alpha_list = c(.5,1,1.5)
K_list = c(3,5,8)

for(p in 2:length(alpha_list)){
  for(pp in 1:length(K_list)){

#Hyperparameters
N=100
a=1
b=1
K_max=4
M=3800
max_clust=3
min_clust= K_list[pp]
max_number_games = 100
N_iter = 5000
alpha= alpha_list[p]

beta_max = .8

iteration_title = paste("K",K_list[pp], "alpha" ,alpha_list[p])
iteration_title= gsub(" ", "",as.character(iteration_title))
iteration_title= gsub(".", "_",as.character(iteration_title),fixed = T)

print(iteration_title)

set.seed(13)

synth_power = simulating_tournament_powerlaw(N=N,alpha = alpha,beta_max = beta_max,min_clust = min_clust,max_clust = max_clust , M = M, n_ij_max =max_number_games )
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
beta_0=1  
sigma0=.2
#proposing a new p
p_current <- p_proposal(rbeta(K_true**2,1,1),2,K = K_true)


p.container[,,1] = p_current



#Setting up quantities needed within computations

# here we are transforming a Nx1 vector, containg labels 1...K into
# NXK matrix. z_mat_current[i,k] =1 if node i is in cluster k
z_mat_current= vec2mat(z_current)

#creating an NxN matrix where p_n_current[i,j] is the probability that player i wins vs player j
matrix_z_p_current = p_current%*%t(z_mat_current)
p_n_current = z_mat_current%*%matrix_z_p_current
p_ij_current =  p_n_current[non_negative_n_ij]


A_seq[1] = sum(dbinom(y_ij, n_ij, p_ij_current, log=T)) + get_B(p_current,beta_0)



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
  for(i in 1:N){
    #proposing a new label
    new_label = sample(x=labels_available,size=1)
    
    #updating labels
    z_prime = z_current
    z_prime[i] = new_label

    #computing the new victory probabilities
    #computing the new victory probabilities
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
    
    matrix_z_prime_p_current = p_current %*% t(z_mat_prime)
    p_n_z_prime = z_mat_prime%*%matrix_z_prime_p_current
    p_ij_z_prime =  p_n_z_prime[non_negative_n_ij]
    
    #acceptance ratio
    r = sum(dbinom(y_ij, n_ij, p_ij_z_prime, log=T)) - sum(dbinom(y_ij, n_ij, p_ij_current, log = T))
    
    alpha = min(1, exp(r))
    u = runif(1)
    
    #if accepted
    if(u<alpha){
      acc.count_z = acc.count_z + 1
      z_current=z_prime
      p_ij_current= p_ij_z_prime
      #if not accepted
    }
    z.container[, j] = z_current
  }
  
  #Update of the P matrix
  #----

  #proposing a new p
  p_prime = p_proposal(p = p_current,sigma_p = sigma0*exp((acc.count_p+1)/j),K = K_true)
  p_prime[is.na(p_prime)] <- 0
  p_prime = (1 -t(p_prime*upper.tri(p_prime))) * (lower.tri(p_prime, diag = F)*1) + upper.tri(p_prime, diag = T) *p_prime
  diag(p_prime)=0.5
  
  
  #Updating p_ij_prime with the last membership
  z_mat_current = vec2mat(z_current)
  matrix_z_p_prime  = p_prime%*%t(z_mat_current)
  p_n_prime = z_mat_current%*%matrix_z_p_prime
  p_ij_prime =  p_n_prime[non_negative_n_ij]
  
  r = (sum(dbinom(y_ij, n_ij, p_ij_prime, log=T)) +  +get_B(p_prime, beta_0)) - 
    (sum(dbinom(y_ij, n_ij, p_ij_current, log = T)) + +get_B(p_current, beta_0) )
  
  alpha = min(1, exp(r))
  u = runif(1)
  if(u<alpha){
    #counting number of accepted proposals
    acc.count_p = acc.count_p+1
    #updating quantities
    p_ij_current = p_ij_prime
    p_current = p_prime
  }
  #storing results for diagnostics
  p.container[,,j] = p_current

  A_seq[j] = sum(dbinom(y_ij, n_ij, p_ij_current, log=T)) +  get_B(p_current,beta_0)
  j=j+1
}


#diagnostic for z
#-----
#proposal performance
acceptance_rate= acc.count_z/(N*N_iter)*100

#mixing
ts.plot(A_seq[-c(1:N_iter*0.5)])
acf(A_seq[-c(1:N_iter*0.5)])


#estimates
similarity_matrix = pr_cc(z.container[,-c(1:N_iter*0.5)])

point_est = minVI(similarity_matrix)$cl
adj.rand.index(point_est, z.true)

# similarity_plot(y_ij_matrix,z.true,z.true)
# similarity_plot(similarity_matrix,z.true,z.true)

#diagnostic for P
#-----

acceptance_rate_p  = (acc.count_p/N_iter)*100



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
# p_combined = patchwork::wrap_plots(plots, ncol = 3, nrow = 3)
# p_combined


png(as.character(paste("adjacencysimple",iteration_title,".png",sep = "")),width = 800, height = 540)
similarity_plot(y_ij_matrix,z.true,z.true)
dev.off()

png(as.character(paste("similaritysimple",iteration_title,".png",sep = "")),width = 800, height = 540)
similarity_plot(similarity_matrix,z.true,z.true)
dev.off()

png(as.character(paste("heatmap_mse_simple",iteration_title,".png",sep = "")),width = 400, height = 400)
heat_map_blue(matrix = mse_table, "heatmap_mse")
dev.off()


mse_iterations[iteration_title] = sqrt(mse_total)
rand_index_list[iteration_title] = adj.rand.index(point_est, z.true)


  }
  }





