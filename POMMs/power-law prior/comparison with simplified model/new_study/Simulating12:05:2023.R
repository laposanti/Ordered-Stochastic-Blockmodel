
library(fossil)
library(dplyr)
#Hyperparameters
N=100
a=1
b=1
K=3

M=8000
max_clust=3
min_clust= 3


N_iter = 10000
alpha=1


gamma_vec = rep(1,min_clust)
beta_max=.8

set.seed(16)
synth = simulating_tournament_powerlaw_new(N,alpha,beta_max,K, M = M,gamma_vec = gamma_vec,n_ij_max = 6)

K_true=K
z.true = synth$z_true
labels_available=1:K

n_ij_matrix = synth$n_ij_true
y_ij_matrix= synth$y_ij_true



#selecting just those values such that there is inside a non-0 entry
upper.tri_n_ij = upper.tri(n_ij_matrix)
non_negative_n_ij = which(upper.tri_n_ij & n_ij_matrix > 0, arr.ind = T)



#retrieving the non-0 values
n_ij = n_ij_matrix[non_negative_n_ij]
y_ij = y_ij_matrix[non_negative_n_ij]


p_ij_true= synth$p_ij_true[non_negative_n_ij]








#Containers for p
#---
alpha.container = matrix(0,N_iter,1)
p.container = array(0,c(K_true,K_true,N_iter))



#Containers for z
#---
z.container = matrix(0,N,N_iter)
A_seq = matrix(0,N_iter,1)


#Containers for p
#---
alpha.container = matrix(0,N_iter,1)
p.container = array(0,c(K_true,K_true,N_iter))


#Initialization for the z vector
#---
init = kmeans(x = y_ij_matrix,centers =K)$cluster
adj.rand.index(init,z.true)
z_current= init
z.container[,1] = z_current

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

#containers for the counts of accepted proposals
acc.count_z = 0
acc.count_p = 0

#setting time tracker
pb=txtProgressBar(min=1,max=N_iter)
j=2

diagnostic = matrix(0, nrow=3,ncol=N_iter)

#READY TO BOMB!
while (j < N_iter + 1) {
  setTxtProgressBar(pb, j)
  
  #Complete sweeep of z vector
  #----
  sweeping_order = sample(x=c(1:N),size=N, replace=F)
  for(i in sweeping_order){
    #proposing a new label
    new_label = adjacent_label_sampler(labels_available = labels_available, z_current[i])
    
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
        new_label = adjacent_label_sampler(labels_available=labels_available, z_current[i])
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
    z.container[, j] = z_current
  }
  
  #Update of the P matrix
  #----
  sigma_prime = 1
  
  #proposing a new alpha
  alpha_prime <- rlnorm_param(1,z.mu = alpha_current, z.sigma =sigma_prime)
  
  #generating a proposal matrix
  p_prime = simulating_POMM_powerlaw(K_true,alpha_prime,beta_max)
  
  
  
  #Updating p_ij_prime with the last membership
  z_mat_current = vec2mat(z_current)
  matrix_z_p_prime  = p_prime$matrix%*%t(z_mat_current)
  p_n_prime = z_mat_current%*%matrix_z_p_prime
  p_ij_prime =  p_n_prime[non_negative_n_ij]
  
  r= (sum(dbinom(y_ij, n_ij, p_ij_prime, log=T))+  l_like_p_ij(p_prime$matrix,p_prime$truncations)+dlnorm_param(alpha_prime))  - 
    (sum(dbinom(y_ij, n_ij, p_ij_current, log = T))+  l_like_p_ij(p_current$matrix,p_current$truncations)+dlnorm_param(alpha_current))
  
  alpha_r = min(1, exp(r))
  u = runif(1)
  if(u<alpha_r){
    diagnostic[1,j-1] = 1
    diagnostic[2,j-1] = alpha_prime
    diagnostic[3,j-1] = mean((p_n_prime- p_n_true**2))
    #counting number of accepted proposals
    acc.count_p = acc.count_p+1
    #updating quantities
    p_ij_current = p_ij_prime
    alpha_current = alpha_prime
    p_current = p_prime
  }else{
    diagnostic[1,j-1] = 0
    diagnostic[2,j-1] = alpha_prime
    diagnostic[3,j-1] = mean((p_n_prime- p_n_true**2))
  }
  #storing results for diagnostics
  p.container[,,j] = p_current$matrix
  alpha.container[j]= alpha_current
  
  A_seq[j] = sum(dbinom(y_ij, n_ij, p_ij_current, log=T)) +  l_like_p_ij(p_current$matrix,p_current$truncations) +
    dlnorm_param(alpha_prime) + ddirichlet_multinomial(N,K_true,n_k = n_k_current,my_alpha = gamma_vec)
  j=j+1
}

#diagnostic for z
#-----
#proposal performance
acceptance_rate= acc.count_z/(N*N_iter)*100
print(acceptance_rate)
#mixing
ts.plot(A_seq[-c(1:N_iter*0.5)], main="Traceplot of p(y_ij|z,P) p(z) p(P)",
        xlab = "Iterations", ylab ="p(y_ij|z,P) p(z) p(P)")

acf(A_seq[-c(1:N_iter*0.5)],main="Autocorrelation plot",
    xlab = "Lag", ylab ="ACF")


#estimates
similarity_matrix = pr_cc(z.container[,-c(1:N_iter*0.25)])



point_est = minVI(similarity_matrix)$cl

adj.rand.index(point_est, z.true)


similarity_plot(y_ij_matrix,z.true,z.true)



similarity_plot(similarity_matrix,z.true,z.true)

#diagnostic for P
#-----

acceptance_rate_p  = (acc.count_p/N_iter)*100
print(acceptance_rate_p)

plot(ts(alpha.container), main = "Traceplot of alpha values")

abline(h = alpha, col = "red", lty = 2)
mean(alpha.container)


for(i in 1:K_true){
  for(j in 1:K_true){
    print(mean(p.container[i,j,]) - synth$P_matrix[i,j])
  }
}

diagnostic_df = data.frame(t(diagnostic))
names(diagnostic_df) = c('y.n',"my_alpha", "err")

ggplot(diagnostic_df, aes(x = my_alpha,y=err, col=factor(y.n)))+
  geom_point()+
  geom_vline(aes(xintercept= alpha), col="red")

ggplot(diagnostic_df, aes(x = my_alpha, col=factor(y.n)))+
  geom_histogram(aes(fill=factor(y.n)))+
  geom_vline(aes(xintercept= alpha), col="red")

ggplot(diagnostic_df, aes(x = log(err), col=factor(y.n)))+
  geom_histogram(aes(fill=factor(y.n)))




burnin_p = p.container[,,-(N_iter*0.5)]

plots = list()
for(i in 1:K_true) {
  for(j in 1:K_true) {
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




# Create a data frame containing the observations
df_alpha <- data.frame(alpha_est = alpha.container)

