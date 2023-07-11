


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
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/functionP_POMM2.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_overlap.R")
N_iter=10000
set.seed(34)

N=100
M= 40000
K=5
alpha=0.5
beta_max= .85
overlap=0.2
gamma_vec = vector()
for(i in 1:K){
  gamma_vec = append(gamma_vec, i/(K**2))
}
model='POMM'
diag0.5<-T

synth = simulating_tournament_new_overlap_norm(N = N, alpha = alpha,overlap = overlap,
                                               beta_max = beta_max,
                                               K=K, M = M,
                                               gamma_vec = gamma_vec,
                                               n_ij_max = 6,model = model,diag0.5 = T 
)

n_ij_matrix=synth$n_ij_true
y_ij_matrix=synth$y_ij_true
p_ij_true=synth$p_ij_true
z_true=synth$z_true
P_true= synth$P_matrix


similarity_plot(y_ij_matrix, z_true,z_true)







barplot(rowSums(y_ij_matrix)/rowSums(n_ij_matrix))
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


similarity_plot(y_ij_matrix,z_true,z_true)

#setting containers
z_container= matrix(0, nrow=N, ncol=N_iter)


#initializing quantities
sigma_prime=.15
alpha_current = 1
truncations_current <- improper_prior5(K,beta_max,alpha = alpha_current,diag0.5 = T)

#generating a proposal matrix
alpha_container = matrix(0, nrow=1, ncol=N_iter)
overlap_container =  matrix(0, nrow=1, ncol=N_iter)
p_container = array(0, dim=c(K,K,N_iter))

p_current = simulating_POMM_powerlaw2(K,alpha_current,truncations_current,beta_max,diag0.5 = T)

overlap_current= runif(1,0.1,0.8)


A_container = matrix(0, nrow=1, ncol=N_iter)
B_container = matrix(0, nrow=1, ncol=N_iter)
C_container = matrix(0, nrow=1, ncol=N_iter)


#initializing_containers
init<-kmeans(x = y_ij_matrix,centers = K)$cluster
z_current=  z_true
adj.rand.index(init,z_true)




n_k_current = as.vector(table(z_current))
z_mat_current = vec2mat(z_current)
#p_ij_function = calculate_victory_probabilities(z_mat_current,P_true)
aux = p_current%*%t(z_mat_current)
p_nbyn_current = z_mat_current%*%aux
p_ij_current = p_nbyn_current[upper.tri.non.zero]


similarity_plot(p_nbyn_current, z_current,z_current)



labels_available = 1:K

A_current= sum(dbinom(y_ij, n_ij, p_ij_current, log = T))
B_current=ddirichlet_multinomial(N,K,n_k = n_k_current ,my_alpha = gamma_vec)
C_current =  l_like_p_ij_normal_overlap(K = K, P_matrix = p_current,overlap =overlap_current, truncations = truncations_current,diag0.5 = T) + + dlnorm_mu_sigma(overlap_current)+ dlnorm_mu_sigma(alpha_current)

z_container[,1] = z_current
p_container[,,1] = p_current
alpha_container[1] = alpha_current
overlap_container[1] = overlap_current


A_container[1]=A_current
B_container[1]=B_current
C_container[1]=C_current
#containers for the counts of accepted proposals
acc.count_z = matrix(1, N, 1)
acc.count_overlap=1
acc.count_p = matrix(1,K,K)
acc.count_alpha = 1

sigma_z <- rep(0.5,N)
sigma_alpha = 0.2
sigma_overlap =0.2 
sigma_p= matrix(rep(0.2, K**2),K,K)

sigma_z_container<- matrix(0, N, N_iter)
sigma_alpha_container <- matrix(0, N_iter,1)
sigma_overlap_container<- matrix(0, N_iter,1)
sigma_p_container<- array(0, dim=c(K,K, N_iter))


#setting time tracker
pb=txtProgressBar(min=1,max=N_iter)
j=2
optimal_p = 0.3





#READY TO BOMB!



for(j in 1:N_iter){
  
  setTxtProgressBar(pb, j)
  # Create a vector of update indices in random order
  updateOrder <- c(1:4)
  
  for (updateIndex in updateOrder) {
    if (updateIndex == 1) {
      #alpha_update
      alpha_update = P_POMM_alpha_update(z_current = z_current,
                                         p_current = p_current,
                                         K = K,n_ij = n_ij,
                                         y_ij = y_ij,
                                         A_current = A_current, C_current = C_current,
                                         upper.tri.non.zero = upper.tri.non.zero,
                                         alpha_current = alpha_current,
                                         beta_max = beta_max,overlap_current = overlap_current,acc.count_alpha = acc.count_alpha, 
                                         sigma_alpha = sigma_alpha, truncations_current = truncations_current,labels_available = labels_available)
      
      #updating quantities
      sigma_alpha = alpha_update$sigma_alpha 
      acc.count_alpha = alpha_update$acc.moves
      C_current = alpha_update$C_current
      alpha_current = alpha_update$alpha_current
      truncations_current = alpha_update$truncations_current
      #adaptive standard deviation
      if(j %% 50 == 0){
        sigma_alpha = tuning_proposal(iteration=j,acceptance_count = acc.count_alpha,sigma = sigma_alpha,acceptanceTarget = optimal_p,min_sigma = 0.02)
      }
    }else if (updateIndex == 2) {
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

      
      #overlap UPDATE----------------------------------------------------------------
      
      }else if (updateIndex == 3) {overlap_update= P_POMM_overlap_update(z_current = z_current,
                                            p_current = p_current,
                                            K = K,n_ij = n_ij,
                                            y_ij = y_ij,
                                            A_current = A_current, C_current = C_current,
                                            upper.tri.non.zero = upper.tri.non.zero,
                                            alpha_current = alpha_current,
                                            beta_max = beta_max,overlap_current = overlap_current,
                                            acc.count_overlap = acc.count_overlap,sigma_overlap = sigma_overlap,truncations_current = truncations_current)
      
      #updating quantities
      sigma_overlap = overlap_update$sigma_overlap
      acc.count_overlap = overlap_update$acc.moves
      C_current = overlap_update$C_current
      overlap_current = overlap_update$overlap_current
      if(j %% 50 == 0){
        sigma_overlap = tuning_proposal(iteration=j,acceptance_count = acc.count_overlap,sigma = sigma_overlap,acceptanceTarget = optimal_p,min_sigma = 0.02)
      }
    }else if (updateIndex == 4) {
      #P UPDATE----------------------------------------------------------------
      
      p_update= P_POMM_update_given_overlap1(z_current = z_current,
                                             p_current = p_current,
                                             K = K,n_ij = n_ij,
                                             y_ij = y_ij,
                                             A_current = A_current, C_current = C_current,
                                             upper.tri.non.zero = upper.tri.non.zero,
                                             alpha_current = alpha_current,
                                             beta_max = beta_max,overlap_current = overlap_current,diag0.5=T,
                                             labels_available = labels_available,
                                             acc.count_p =acc.count_p,sigma_p = sigma_p,truncations_current = truncations_current)
      
      #updating quantities
      sigma_p = p_update$sigma_p
      C_current = p_update$C_current
      A_current = p_update$A_current
      p_current = p_update$p_current
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
    }
  }
    if(j %% 1000 == 0){
      print( paste("iteration",j,"acceptance", acc.count_p))
    }
  
    #storing scales
    sigma_alpha_container[j] <- sigma_alpha
    sigma_overlap_container[j]<- sigma_overlap
    sigma_p_container[,,j]<- sigma_p
    sigma_z_container[,j] <- sigma_z
    
    #storing results for inference
    A_container[j] = A_current
    B_container[j] = B_current
    C_container[j]= C_current
    
    alpha_container[j]<- alpha_current
    z_container[,j] <- z_current
    overlap_container[j] = overlap_current
    p_container[,,j] = p_current
  
}



MAP_retriever = A_container + B_container + C_container
z_MAP= z_container[,which(A_container == max(A_container))]

similarity_matrix = pr_cc(z_container[,-c(1:N_iter*0.55)])
point_est = minVI(similarity_matrix)$cl
adj.rand.index(point_est, z_true)


similarity_plot(y_ij_matrix,synth$z_true,synth$z_true)
similarity_plot(similarity_matrix,synth$z_true,synth$z_true)

burnin_p = p_container[,,-(N_iter*0.5)]





plot(ts(sigma_alpha_container))


burnin_p = p_container[,,-(N_iter*0.5)]

ts.plot(A_container[-c(1:5000)])
acf((A_container[-c(1:5000)]))

ts.plot(p_container[1,3,-c(1:100)])


ts.plot(C_container[-c(1:8000)])
acf((C_container[-c(1:8000)]))



ts.plot(t(overlap_container))
mean(overlap_container[-c(1:(N_iter*0.5))])
example=data.frame(x1 = t(overlap_container))
ggplot(example, aes(x= x1))+
  geom_density()+
  geom_vline(xintercept = overlap)

burnin_p = p_container[,,-c(1:(N_iter*0.5))]

ts.plot(t(alpha_container))
mean(alpha_container[-c(1:(N_iter*0.5))])









alpha_container[which(C_container == max(C_container))]


burnin_p = p_container[,,-(N_iter*0.25)]
#acceptance rate
acc.count_p/( N_iter)


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

example = data.frame(y=(burnin_p[1,3,]))
ggplot(example, aes(x=y))+
  geom_density()+
  geom_vline(xintercept =synth$P_matrix[1,3])

pander::pander(mse_table)


P_est = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    P_est[i,j]= mean(burnin_p[i,j,]) 
  }
}
pander::pander(P_est)

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
for (i in 1:N_new){
  for (j in 1:N){
    Yij_new[i,j] <-rbinom(1,1,prob= synth$P_matrix[z_new[i],z_true[j]])
  }
}



z_proportion = table(point_est)/N


z_predicted_prob = matrix(0, N_new, K)
for(k in labels_available){
  for(i in 1:N_new){
    lik_i  <- 0 
    for(j in 1:N){
      lik_i  <- lik_i + dbinom(Yij_new[i,j],1, P_est[k, point_est[j]],log = T)
    }
    lik_i <- lik_i + log(z_proportion[k])
    z_predicted_prob[i,k] = lik_i 
  }
}



misclassification_rate <- (N_new-sum(diag(table(z_new,apply(z_predicted_prob,1,which.max)))))/N_new

print(misclassification_rate)

