library(foreach)
library(truncnorm)
library(doParallel)
library(doSNOW)

source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/function_Z1.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/Modular_code/functionP_POMM2.R")
N_iter=20000
set.seed(34)

N=100
M= 10000
K=4
alpha=1

beta_max= .85



gamma_vec = vector()
for(i in 1:K){
  gamma_vec = append(gamma_vec, i/(K**2))
}


synth = simulating_tournament_new_norm(N = N, alpha = alpha,
                                       beta_max = beta_max,
                                       K=K, M = M,
                                       gamma_vec = gamma_vec,
                                       n_ij_max = 6,model = 'POMM',diag0.5 = T 
)

n_ij_matrix=synth$n_ij_true
y_ij_matrix=synth$y_ij_true
p_ij_true=synth$p_ij_true
z_true=synth$z_true
P_true= synth$P_matrix






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



# Set the number of cores to use
num_cores <- 4  # Adjust this number based on your available cores

initial_z<- matrix(0,N,num_cores)
initial_p<- array(0,dim=c(K,K,num_cores))
initial_A<- matrix(0,1,num_cores)
initial_B<- matrix(0,1,num_cores)
initial_C<- matrix(0,1,num_cores)
initial_alpha<- matrix(0,1,num_cores)


for(i in 1:num_cores){
  #initializing quantities
  sigma_prime=.15
  alpha_current = 1
  truncations_current <- improper_prior5(K,beta_max,alpha = alpha_current)
  
  #generating a proposal matrix
  alpha_container = matrix(0, nrow=1, ncol=N_iter)
  p_container = array(0, dim=c(K,K,N_iter))
  p_current = simulating_POMM_powerlaw2(K,alpha_current,truncations_current,beta_max)
  initial_p[,,i]<-p_current
  


#initializing quantities
sigma_prime=.15
alpha_current = runif(1,0.5,1.5)
initial_alpha[i]<-alpha_current
truncations_current <- improper_prior5(K,beta_max,alpha = alpha_current)

#generating a proposal matrix
alpha_container = matrix(0, nrow=1, ncol=N_iter)
p_container = array(0, dim=c(K,K,N_iter))
p_current = simulating_POMM_powerlaw2(K,alpha_current,truncations_current,beta_max)


A_container = matrix(0, nrow=1, ncol=N_iter)
B_container = matrix(0, nrow=1, ncol=N_iter)
C_container = matrix(0, nrow=1, ncol=N_iter)
acc.count_p = 0

#initializing_containers
init<-kmeans(x = y_ij_matrix,centers = K)$cluster
z_current= init
adj.rand.index(init,z_true)

initial_z[,i]<- init


n_k_current = as.vector(table(z_current))
z_mat_current = vec2mat(z_current)
#p_ij_function = calculate_victory_probabilities(z_mat_current,P_true)
aux = p_current%*%t(z_mat_current)
p_nbyn_current = z_mat_current%*%aux
p_ij_current = p_nbyn_current[upper.tri.non.zero]




labels_available = 1:K

A_current= sum(dbinom(y_ij, n_ij, p_ij_current, log = T))
B_current=ddirichlet_multinomial(N,K,n_k = n_k_current ,my_alpha = gamma_vec)
C_current =  l_like_p_ij_normal(K = K, P_matrix = p_current,truncations = truncations_current,diag0.5 = T) + dlnorm_param(alpha_current)

initial_A[i]<- A_current
initial_B[i]<- B_current
initial_C[i]<- C_current


}





# Set the number of iterations
N_iter <-20000

cl <- makeCluster(num_cores)
registerDoParallel(cl)
# Create a list to store results for each chain
results <- foreach(chain_index = 1:num_cores, .packages = c("progress")) %dopar% {
  library(truncnorm)
  # Initialize the variables for this chain
  z_current <- initial_z[,chain_index]  # Replace initial_z with your initial values for z
  p_current <- initial_p[,,chain_index]  # Replace initial_p with your initial values for p
  A_current <- initial_A[chain_index]  # Replace initial_A with your initial values for A
  B_current <- initial_B[chain_index]  # Replace initial_B with your initial values for B
  C_current <- initial_C[chain_index]  # Replace initial_C with your initial values for C
  alpha_current <- initial_alpha[chain_index]  # Replace initial_alpha with your initial values for alpha
  
  acc.count_z <- 0
  acc.count_p <- 0
  
  z_container <- matrix(0, nrow = N, ncol = N_iter)
  A_container <- numeric(N_iter)
  B_container <- numeric(N_iter)
  C_container <- numeric(N_iter)
  alpha_container <- numeric(N_iter)
  p_container <- array(0, dim = c(K, K, N_iter))
  
  #pb <- pbapply::pbapply(1:N_iter, cl = cl, label = paste("Chain", chain_index, "Progress"), level = 1, silent = FALSE)
  
  for (j in 1:N_iter) {
    # Z UPDATE ----------------------------------------------------------------
    z_sweep <- z_update_1(z_current, A_current, B_current, y_ij, n_ij, p_current, labels_available = labels_available, upper.tri.non.zero = upper.tri.non.zero, gamma_vec = gamma_vec, K = K)
    
    acc.count_z <- acc.count_z + z_sweep$acc.moves
    z_current <- z_sweep$z_current
    n_k_current <- z_sweep$n_k_current
    B_current <- z_sweep$B_current
    A_current <- z_sweep$A_current
    
    z_container[, j] <- z_current
    
    # P UPDATE ----------------------------------------------------------------
    p_update <- P_POMM_update2(z_current = z_current,
                               p_current = p_current,
                               K = K, n_ij = n_ij, y_ij = y_ij,
                               A_current = A_current, C_current = C_current,
                               upper.tri.non.zero = upper.tri.non.zero,
                               alpha_current = alpha_current,
                               beta_max = beta_max)
    
    acc.count_p <- acc.count_p + p_update$acc.moves
    
    p_current <- p_update$p_current
    C_current <- p_update$C_current
    A_current <- p_update$A_current
    alpha_current <- p_update$alpha_current
    
    # Store results for inference
    A_container[j] <- A_current
    B_container[j] <- B_current
    C_container[j] <- C_current
    alpha_container[j] <- alpha_current
    p_container[, , j] <- p_current
    
    #pb$tick()  # Update the progress bar
  }
  
  # Return the results for this chain
  list(z_container = z_container, A_container = A_container, B_container = B_container, C_container = C_container, alpha_container = alpha_container, p_container = p_container)
  
  }

# Stop the parallel backend
stopCluster(cl)
registerDoSEQ()

# Access the results for each chain
chain_1_results <- results[[1]]
chain_2_results <- results[[2]]
chain_3_results <- results[[3]]
chain_4_results <- results[[4]]





#Diagnostics for the P entries ---------
entry1=1
entry2=2
#checking convergence
density_plot_df<-data.frame(chain1_density = chain_1_results$p_container[entry1,entry2,],
           chain2_density = chain_2_results$p_container[entry1,entry2,],
           chain3_density = chain_3_results$p_container[entry1,entry2,],
           chain4_density = chain_4_results$p_container[entry1,entry2,])


# Define the blue shades for the color gradient
blue_shades <- c("#77B7D1", "#3689A8", "#22677F", "#134B5C")

# Create a density plot of points for each level
ggplot() +
  # Add a layer for each level
  lapply(seq(1,4), function(i) {
    geom_density(data = data.frame(x = density_plot_df[,i]), aes(x = x, y = ..density.., color = paste("Chain ", i)), alpha = .5)
  }) +
  geom_vline(aes(xintercept = P_true[entry1,entry2]), col="red", alpha=0.5)+
  # Set the x-axis limits
  scale_x_continuous(limits = c(min(density_plot_df), max(density_plot_df))) +# Remove the legend title
  guides(color = guide_legend(title = NULL)) +
  # Set the legend title
  labs(fill = "Chain", x = "Points", y = "Density", title = paste("Density Plot of Entries","P[",entry1,",",entry2,"]"))+
  # Set blue color scale for the lines
  scale_color_manual(values=blue_shades)+
  theme_bw()








#Traceplots likelihood -------------------

chain_1_A <- chain_1_results$A_container
chain_2_A <- chain_2_results$A_container
chain_3_A <- chain_3_results$A_container
chain_4_A <- chain_4_results$A_container


ts.plot(chain_1_A[-(1:12000)], main= "Traceplot Likelihood Chain 1")
ts.plot(chain_2_A[-(1:12000)], main= "Traceplot LikelihoodChain 2")
ts.plot(chain_3_A[-(1:12000)], main= "Traceplot LikelihoodChain 3")
ts.plot(chain_4_A[-(1:12000)], main= "Traceplot LikelihoodChain 4")



#autocorrplot

acf(chain_1_A[-(1:12000)], main= "Autocorrplot Likelihood Chain 1")
acf(chain_2_A[-(1:12000)], main= "Autocorrplot LikelihoodChain 2")
acf(chain_3_A[-(1:12000)], main= "Autocorrplot LikelihoodChain 3")
acf(chain_4_A[-(1:12000)], main= "Autocorrplot LikelihoodChain 4")



#posterior predictive checks

z_1_container <- chain_1_results$z_container[,-c(1:N_iter*0.5)]
z_2_container <- chain_2_results$z_container[,-c(1:N_iter*0.5)]
z_3_container <- chain_3_results$z_container[,-c(1:N_iter*0.5)]
z_4_container <- chain_4_results$z_container[,-c(1:N_iter*0.5)]

P_k_container_1 <- chain_1_results$p_container[,,-c(1:N_iter*0.5)]
P_k_container_2 <- chain_2_results$p_container[,,-c(1:N_iter*0.5)]
P_k_container_3 <- chain_3_results$p_container[,,-c(1:N_iter*0.5)]
P_k_container_4 <- chain_4_results$p_container[,,-c(1:N_iter*0.5)]

similarity_matrix = pr_cc(z_4_container)
point_est = minVI(similarity_matrix)$cl

P_zizj = array(0, dim=c(N,N,N_iter*0.5))

for(i in 1:(N_iter*0.25)){
  z_mat_ith = vec2mat(z_4_container[,i])
  P_zizj[,,i] = calculate_victory_probabilities(z_mat_ith,P_k_container_4[,,i])
}

y_ij_sim = array(0, dim=c(N,N,N_iter*0.25))
for(t in 1:(N_iter*0.25)){
  for(i in 1:N){
    for(j in 1:N){
      y_ij_sim[i,j,t] = rbinom(1, n_ij_matrix[i,j], P_zizj[i,j,t])
    }
  }
}

y_ij_sum = matrix(0, N, N_iter*0.25)
for(t in 1:N_iter*0.25){
  y_ij_sum[,t]<- rowSums(y_ij_sim[,,t])
}

df_plot = data.frame(id_players=c(1:N),
                     actual_victories = rowSums(y_ij_matrix), 
                     mean_est_victories = apply(y_ij_sum,1,mean),
                     quantile_0.05 = apply(y_ij_sum,1,quantile,0.05),
                     quantile_0.95 = apply(y_ij_sum,1,quantile,0.95))

blue_shades <- c("#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f")  # Blue color shades for legend


ggplot(df_plot, aes(x = id_players, y = actual_victories, color="label")) +
  geom_point(aes(color = "label"), alpha = 0.8) +
  geom_errorbar(aes(ymin = quantile_0.05, ymax = quantile_0.95, color = "Error Bar"), width = 0.5, alpha = 0.5) +
  geom_point(aes(y = mean_est_victories, color = "Mean Estimation"), size = 1) +
  scale_x_discrete() +
  scale_color_manual(values = blue_shades, labels = c("Error Bar", "Mean Estimation", "Actual Victories")) +
  labs(x = "Players", y = "Number of Victories", title = "Posterior Predictive Check") +
  theme_bw()





