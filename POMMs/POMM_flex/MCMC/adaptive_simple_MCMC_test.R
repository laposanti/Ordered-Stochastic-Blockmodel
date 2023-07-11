
library(mcclust.ext)
library(loo)
library(dbscan)
library(randnet)
library(fossil)
library(dplyr)
library(truncnorm)
library(EnvStats)

source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/function_Z1.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/function_P1.R")


N_iter=10000
set.seed(79141)

N=100
M= 4000
K=5
alpha=0.5

beta_max= .85



gamma_vec = vector()
for(i in 1:K){
  gamma_vec = append(gamma_vec, i/(K**2))
}

diag0.5<-T
synth = simulating_tournament_new_overlap_norm(N = N, alpha = alpha,
                                               beta_max = beta_max,
                                               K=K, M = M,
                                               gamma_vec = gamma_vec,
                                               n_ij_max = 6,model = 'Simple',diag0.5 = diag0.5, overlap = overlap
)


n_ij_matrix=synth$n_ij_true
y_ij_matrix=synth$y_ij_true
p_ij_true=synth$p_ij_true
z_true=synth$z_true
P_true= synth$P_matrix

similarity_plot(y_ij_matrix, z_true,z_true)

barplot(rowSums(n_ij_matrix))
barplot(colSums(y_ij_matrix)/colSums(n_ij_matrix))



#upper.tri.non.zero= which(n_ij_matrix >= 0)
#activate for simple model
upper.tri.non.zero = which(n_ij_matrix > 0,arr.ind = T)
#working with just the upper triangular entries and with non-zero-values
n_ij = n_ij_matrix[upper.tri.non.zero]
y_ij = y_ij_matrix[upper.tri.non.zero]
p_ij = p_ij_true[upper.tri.non.zero]
# n_ij = n_ij_matrix
# y_ij = y_ij_matrix
# p_ij = p_ij_true



#------
# setting containers

z_container= matrix(0, nrow=N, ncol=N_iter)
p_container = array(0, dim=c(K,K,N_iter))
A_container = matrix(0, nrow=1, ncol=N_iter)
B_container = matrix(0, nrow=1, ncol=N_iter)
C_container = matrix(0, nrow=1, ncol=N_iter)

#--------
#initializing quantities

p_current= matrix(runif(K**2,1-beta_max,beta_max),K,K)
diag(p_current) <- rep(0.5,K)
p_current = semi_symmetric(p_current)

z_current= z_true
n_k_current = as.vector(table(z_current))
z_mat_current = vec2mat(z_current)

aux = p_current%*%t(z_mat_current)
p_nbyn_current = z_mat_current%*%aux
p_ij_current = p_nbyn_current[upper.tri.non.zero]

A_current= sum(dbinom(y_ij, n_ij, p_ij_current, log = T))
B_current=ddirichlet_multinomial(N,K,n_k = n_k_current ,my_alpha = gamma_vec)
C_current =  sum(dunif(p_current[upper.tri(p_current)],min = 1-beta_max, max = beta_max, log = T))



labels_available = 1:K

#---------
acc.count_z = rep(1,N)
acc.count_p = matrix(1,K,K)
#updating containers
z_container[,1] = z_current
p_container[,,1] = p_current
A_container[1]=A_current
B_container[1]=B_current
C_container[1]=C_current

#containers for the counts of accepted proposals
sigma_z <- rep(0.5,N)
sigma_p= matrix(rep(0.2, K**2),K,K)

sigma_z_container<- matrix(0, N, N_iter)
sigma_p_container<- array(0, dim=c(K,K, N_iter))

#-----------

#setting time tracker
pb=txtProgressBar(min=1,max=N_iter)
j=2

optimal_p =0.25
#READY TO BOMB!
while (j < N_iter + 1) {
  setTxtProgressBar(pb, j)
  #z UPDATE----------------------------------------------------------------
  
  z_update = z_update_adaptive( z_current = z_current,
                                P_matrix = p_current,
                                K = K,n_ij = n_ij,
                                y_ij = y_ij,
                                A_current = A_current,B_current=B_current,
                                upper.tri.non.zero = upper.tri.non.zero,labels_available = labels_available,
                                gamma_vec = gamma_vec,acc.count_z = acc.count_z,sigma_z = sigma_z)
  
  
  acc.count_z = z_update$acc.moves
  if(j %% 50 == 0){
    for(t in 1:N){
      sigma_z[t] = tuning_proposal(iteration=j,acceptance_count = acc.count_z[t],sigma = sigma_z[t],acceptanceTarget = optimal_p,min_sigma = 0.2)
    }
    
  }
  #updating quantities
  z_current <- z_update$z_current
  B_current<- z_update$B_current
  A_current = z_update$A_current
  
  #A_seq[j]= sum(dbinom(y_ij, n_ij, p_ij_current, log = T)) + ddirichlet_multinomial(N,K_true,n_k = n_k_current ,my_alpha = gamma_vec)
  z_container[, j] = z_current
  
  #P UPDATE----------------------------------------------------------------
  
  p_update= P_simple_update_adaptive(z_current = z_current,
                                     p_current  = p_current,K = K,n_ij = n_ij,y_ij = y_ij,diag0.5 = diag0.5,
                                     A_current = A_current,C_current = C_current,
                                     upper.tri.non.zero = upper.tri.non.zero,acc.count_p =  acc.count_p,sigma_p = sigma_p,beta_max=beta_max,labels_available = labels_available)
  
  acc.count_p = p_update$acc.moves
  if(j %% 50 == 0){
    j_start = ifelse(diag0.5, yes = 1, no = 0)
    K_stop = ifelse(diag0.5, yes = K-1, no = K)
    for( ii in 1:K_stop){
      for(jj in (ii+j_start):K){
        sigma_p[ii,jj] <- tuning_proposal(iteration=j,acceptance_count = acc.count_p[ii,jj],sigma = sigma_p[ii,jj],acceptanceTarget = optimal_p,min_sigma = 0.005)
      }
    }
  }
  
  #storing scales
  sigma_p_container[,,j]<- sigma_p
  sigma_z_container[,j] <- sigma_z
  
  #updating quantities
  p_current=  p_update$p_current
  A_current = p_update$A_current
  C_current = p_update$C_current
  
  #storing results for inference
  A_container[j] = A_current
  B_container[j] = B_current
  C_container[j]= C_current
  p_container[,,j] = p_current
  
  j=j+1
}



ts.plot(A_container[-c(1:20)])
acf((A_container[-c(1:10)]))


ts.plot(B_container[-c(1:100)])
acf((B_container[-c(1:100)]))


posterior_prop = A_container+B_container

ts.plot(posterior_prop[-c(1:100)])
MAP = z_container[, which(posterior_prop == max(posterior_prop))[1]]



acceptance_rate= acc.count_z/(N*N_iter)*100
print(acceptance_rate)

#mixing
# ts.plot(A_seq[-c(1:N_iter*0.5)], main="Traceplot of p(y_ij|z,P) p(z) p(P)",
#         xlab = "Iterations", ylab ="p(y_ij|z,P) p(z) p(P)")
# dev.off()
# acf(A_seq[-c(1:N_iter*0.5)],main="Autocorrelation plot",
#     xlab = "Lag", ylab ="ACF")
# dev.off()


#estimates
similarity_matrix = pr_cc(z_container[,-c(1:N_iter*0.50)])
point_est = minVI(similarity_matrix)$cl

adj.rand.index(MAP, z_true)
adj.rand.index(point_est, z_true)



similarity_plot(y_ij_matrix,synth$z_true,synth$z_true)
similarity_plot(similarity_matrix,synth$z_true,synth$z_true)


burnin_p = p_container[,,-(N_iter*0.25)]
burnin_posterior = A_container[-c(1:5000)]

ts.plot(burnin_posterior)
summary(burnin_posterior)





#acceptance rate
acc.count_p/(K**2 * N_iter)


plots = list()
for(i in 1:K) {
  for(j in 1:K) {
    y_try = data.frame(y = as.vector(burnin_p[i, j,]))
    p1 = ggplot(y_try, aes(y)) +
      geom_density(fill = "dodgerblue", alpha = 0.5) +
      scale_x_log10() +
      geom_vline(xintercept = synth$P_matrix[i, j], color = "red")+
      xlab("probability") +
      ylab("Density") +
      ggtitle(paste("Density plot of entry ", i, ",", j, sep = ""))
    
    plots[[length(plots) + 1]] <- p1
  }
}
p_combined = patchwork::wrap_plots(plots, ncol = K, nrow = K)
p_combined

mse_table = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    mse_table[i,j]= (mean(burnin_p[i,j,]) - synth$P_matrix[i,j])
  }
}
pander::pander(mse_table)


P_est = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    P_est[i,j]= mean(burnin_p[i,j,]) 
  }
}
pander::pander(P_est)

vi.dist(point_est,z_true)

vi.dist(MAP,z_true)

# ------------------------------------
# LOUVAIN ALGORITHM
# ------------------------------------
library(igraph)
# transform the adjacency matrix into an igraph object
net <- graph.adjacency(y_ij_matrix, mode=c("directed"),diag = FALSE)

# point estimate
Louv <- cluster_edge_betweenness(net)
# estimated H
length(table(Louv$membership))
# VI distance between estimated and true partition
vi.dist(z_true,t(Louv$membership))

# ------------------------------------
# Regularised SPECTRAL CLUSTERING
# ------------------------------------

reg_sp <- randnet::reg.SP(A = y_ij_matrix,K = 3,iter.max = 1000000000)$cluster
mcclust::vi.dist(z_true,reg_sp)
adj.rand.index(reg_sp, z_true)


# ------------------------------------
# DBSCAN ALGORITHM
# ------------------------------------



# ------------------------------------
# K-MEANS ALGORITHM
# ------------------------------------

# ------------------------------------
# ADJ Rand Index estimate
# ------------------------------------

adj.rand.index(kmeans(y_ij_matrix,3)$cluster,z_true) 
vi.dist(point_est,z_true)

# ------------------------------------
# WAIC estimate
# ------------------------------------

upper.tri.non.zero.waic = which(upper.tri(n_ij_matrix) & n_ij_matrix>0, arr.ind = T)
waic_matrix_container  = matrix(0, nrow=N_iter, ncol=length(upper.tri.non.zero.waic))

for(ii in 1:N_iter){
  z_mat_waic = vec2mat(z_container[,ii])
  p_ij_waic= calculate_victory_probabilities(z_mat_waic, p_container[,,ii])
  waic_matrix_container[ii,]= dbinom(y_ij_matrix[upper.tri.non.zero.waic], size = n_ij_matrix[upper.tri.non.zero.waic], p_ij_waic[upper.tri.non.zero.waic])
}
#measures
waic.matrix(waic_matrix_container)






# --------------------------------------------
# New incoming nodes: Misclassification Error
# --------------------------------------------

N_new = 100 
z_new = c(rep(1,30),rep(2,30),rep(3,40))

# create empty matrix of edges between the 300 new nodes and those in the original network
Yij_new <- matrix(0,N_new,N)



# simulate the new edges
# for (i in 1:N_new){
#   for (j in 1:N){
#     Yij_new[i,j] <-rbinom(1,1,prob= synth$P_matrix[z_new[i],z_true[j]])
#   }
# }

sampled_games = sample(1:N, 100, F)

for (i in 1:N_new){
  for (j in sampled_games){
    Yij_new[i,j] <-rbinom(1,1,prob= synth$P_matrix[z_new[i],z_true[j]])
  }
}


z_proportion = table(point_est)/N


z_predicted_prob = matrix(0, N_new, K)
for(k in labels_available){
  for(i in 1:N_new){
    lik_i  <- 0 
    for(j in sampled_games){
      lik_i  <- lik_i + dbinom(Yij_new[i,j],1, P_est[k, z_MAP[j]],log = T)
    }
    lik_i <- lik_i + log(z_proportion[k])
    z_predicted_prob[i,k] = lik_i 
  }
}



misclassification_rate <- (N_new-sum(diag(table(z_new,apply(z_predicted_prob,1,which.max)))))/N_new

print(misclassification_rate)

calculate_misclassification_rate <- function(N_new, z_new, N, Yij_new, p_est, z_true, sampled_games, labels_available, P_est, z_MAP) {
  # create empty matrix of edges between the N_new nodes and those in the original network
  Yij_new <- matrix(0, N_new, N)
  
  # simulate the new edges
  for (i in 1:N_new){
    for (j in sampled_games){
      Yij_new[i, j] <- rbinom(1, 1, prob = synth$P_matrix[z_new[i], z_true[j]])
    }
  }
  
  z_proportion <- table(point_est) / N
  
  K <- length(labels_available)
  z_predicted_prob <- matrix(0, N_new, K)
  
  for (k in labels_available){
    for (i in 1:N_new){
      lik_i <- 0 
      for (j in sampled_games){
        lik_i <- lik_i + dbinom(Yij_new[i, j], 1, P_est[k, z_MAP[j]], log = TRUE)
      }
      lik_i <- lik_i + log(z_proportion[k])
      z_predicted_prob[i, k] <- lik_i 
    }
  }
  
  misclassification_rate <- (N_new - sum(diag(table(z_new, apply(z_predicted_prob, 1, which.max))))) / N_new
  
  return(misclassification_rate)
}

