
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


N_iter = 20000
alpha=1


gamma_vec = c(1/3, 1/3 ,1/3)
beta_max=.8




set.seed(20)
synth = simulating_tournament_powerlaw_new(N = N,
                                           alpha = alpha,
                                           beta_max = beta_max,
                                           K=K, M = M,
                                           gamma_vec = gamma_vec,
                                           n_ij_max = 6)

K_true=K
z.true = synth$z_true
labels_available=1:K


n_ij_matrix = synth$n_ij_true
y_ij_matrix= synth$y_ij_true



#selecting just those values such that there is inside a non-0 entry
upper.tri_n_ij = upper.tri(n_ij_matrix)
non_negative_n_ij = which((n_ij_matrix > 0) & upper.tri_n_ij, arr.ind = T)



#retrieving the non-0 values
# n_ij = n_ij_matrix[p_logical_mat]
# y_ij = y_ij_matrix[p_logical_mat]
# 
# p_ij_true= synth$p_ij_true[p_logical_mat]
#p_ij_true= synth$p_ij_true[non_negative_n_ij]

# 
# cor(p_ij_true ,y_ij)
# table(p_ij_true[non_negative_n_ij]== p_ij_true[p_logical_mat])
# 
# 
# p_logical_mat = (upper.tri(p_n_current,diag = T) & n_ij_matrix>0)
# maty.n=p_logical_mat
# # Convert matrix to long format data.frame
# df <- reshape2::melt(maty.n)
# names(df) <- c("var1", "var2", "y.n")
# 
# # Create heatmap
# ggplot(df, aes(x = var2, y = desc(var1), fill = factor(y.n))) +
#   geom_tile(color = "white") +
#   scale_fill_manual(values = c("white", "red")) +
#   labs(x = "", y = "") +
#   theme_bw()
# 
# 
# df_1 = data.frame(cbind(non_negative_n_ij, rep(T,3400)))
# names(df_1) <- c("var1", "var2", "y.n")
# # Create heatmap
# ggplot(df_1, aes(x = var2, y = desc(var1), fill = factor(y.n))) +
#   geom_tile(color = "white") +
#   scale_fill_manual(values = c("red")) +
#   labs(x = "", y = "") +
#   theme_bw()

# Merge the two data frames by var1 and var2
# df_merged <- left_join(df, df_1, by = c("var1", "var2"))
# df_merged <- df_merged %>% mutate(y.n.y = ifelse(is.na(y.n.y), 0, 1))  # Replace NAs with 0
# 
# # Create new column to indicate color
# df_merged$color <- ifelse(df_merged$y.n.x == 0 & df_merged$y.n.y == 1 | df_merged$y.n.x == 1 & df_merged$y.n.y == 0 ,"error",
#                           ifelse((df_merged$y.n.x + df_merged$y.n.y== 2), "1", "0"))
# 
# 
# 
# # Create heatmap
# ggplot(df_merged, aes(x = var2, y = desc(var1), fill = color)) +
#   geom_tile(color = "white") +
#   scale_fill_manual(values = c("green", "red","black")) +
#   labs(x = "", y = "") +
#   theme_bw()

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
alpha_current = runif(1,0.01,2)
p_current = simulating_POMM_powerlaw(K_true,alpha_current,beta_max = beta_max)
truncations_current = p_current$truncations
alpha.container[1] = alpha_current
p.container[,,1] = p_current$matrix

#Setting up quantities needed within computations

# here we are transforming a Nx1 vector, containg labels 1...K into
# NXK matrix. z_mat_current[i,k] =1 if node i is in cluster k
z_mat_current= vec2mat(z.true)
n_k_current = colSums(z_mat_current)

#creating an NxN matrix where p_n_current[i,j] is the probability that player i wins vs player j
matrix_z_p_current = p_current$matrix%*%t(z_mat_current)
p_n_current = z_mat_current%*%matrix_z_p_current
p_ij_current = p_n_current[non_negative_n_ij]

if(l_like_p_ij(p_current$matrix,truncations = truncations_current) == -Inf)
{print("ERROR!!")}else{print("p_current - OK")}

matrix_z_p_true = synth$P_matrix%*%t(z_mat_current)
p_n_true = z_mat_current%*%matrix_z_p_true




#p_ij_current =  p_n_current[non_negative_n_ij]

#containers for the counts of accepted proposals
acc.count_z = 0
acc.count_p = 0

#setting time tracker
pb=txtProgressBar(min=1,max=N_iter)
j=2

diagnostic = matrix(0, nrow=5,ncol=N_iter)
entries = expand.grid(c(1:K_true), c(1:K_true))

my_upper_tri = upper.tri(p_current$matrix)
my_lower_tri = lower.tri(p_current$matrix)
#READY TO BOMB!
while (j < N_iter + 1) {
  setTxtProgressBar(pb, j)
  
 
  
  #proposing a new alpha
  alpha_prime <- sample_norm_trunc(1, alpha_current,s =sigma_prime,a = 0.01,b = 3)
  truncations_prime <- improper_prior5(K,beta_max,alpha = alpha_prime)
  #Update of the P matrix
  #----
  sigma_prime = .5
  #generating a proposal matrix
  p_prime = simulating_POMM_powerlaw1(K_true,alpha_prime,truncations_prime,beta_max)
  
  sampled_order = sample(c(0:(K-1)),K,replace = F)
  for(iii in sampled_order){    
    p_scanning = p_current$matrix
    #recovering entries for diag iii for iii=0,...,K_true-1
    entries_diag = cbind(entries[which(entries$Var2 -entries$Var1 == iii),])
  
    #substituing just the diagonal iii
    
    p_scanning[cbind(entries_diag[,1],entries_diag[,2])] = p_prime$matrix[cbind(entries_diag[,1],entries_diag[,2])]
    p_scanning[my_lower_tri] = 1 - t(p_scanning[my_upper_tri])
    
    #tranforming data
    matrix_z_p_scanning  = p_scanning%*%t(z_mat_current)
    p_n_scanning = z_mat_current%*%matrix_z_p_scanning
    p_ij_scanning =  p_n_scanning[non_negative_n_ij]
    
    truncations_scanning = truncations_current
    truncations_scanning[c(iii+1,iii+2)] = truncations_prime[c(iii+1,iii+2)]
    
    r= (sum(dbinom(y_ij, n_ij, p_ij_scanning, log=T))+   l_like_p_ij(p_scanning,truncations_scanning) + dlnorm_param(alpha_prime))  -
    (sum(dbinom(y_ij, n_ij, p_ij_current, log = T))+  l_like_p_ij(p_current$matrix,truncations_current) + dlnorm_param(alpha_current))

  
    alpha_r = min(1, exp(r))
    u = runif(1)
    if(r == -Inf)
    {}else if(u<alpha_r){
      diagnostic[1,j-1] = 1
      diagnostic[2,j-1] = alpha_prime
      diagnostic[3,j-1] = mean((p_n_scanning- p_n_true)**2)
      diagnostic[4,j-1] = mean((p_ij_scanning- p_ij_true)**2)
      diagnostic[5,j-1] = exp(r)
      
      #counting number of accepted proposals
      acc.count_p = acc.count_p+1
      #updating quantities

      alpha_current = alpha_prime
     
      p_current$matrix[cbind(entries_diag[,1],entries_diag[,2])]=p_prime$matrix[cbind(entries_diag[,1],entries_diag[,2])]
      p_current$matrix[my_lower_tri] = 1 - t(p_current$matrix[my_upper_tri])
      #updating truncations
      truncations_current  = truncations_scanning 
      #updating also p_ij
      matrix_z_p_current = p_current$matrix%*%t(z_mat_current)
      p_n_current = z_mat_current%*%matrix_z_p_current
      p_ij_current = p_n_current[non_negative_n_ij]
    }else{
      diagnostic[2,j-1] = alpha_prime
      diagnostic[3,j-1] = mean((p_n_scanning- p_n_true)**2)
      diagnostic[4,j-1] = mean((p_ij_scanning- p_ij_true)**2)
      diagnostic[5,j-1] = exp(r)
    }
  }
  
  #storing results for diagnostics
  p.container[,,j] = p_current$matrix
  alpha.container[j]= alpha_current
  
  A_seq[j] = sum(dbinom(y_ij, n_ij, p_ij_current, log=T)) +  l_like_p_ij(p_current$matrix,truncations = truncations_current) +
    dlnorm_param(alpha_current) + ddirichlet_multinomial(N,K_true,n_k = n_k_current,my_alpha = gamma_vec)
  j=j+1
}

A_burn = A_seq[-c(1:10000)]
alpha_burn = alpha.container[-c(1:10000)]

ts.plot(A_burn,main = "Traceplot of join probability distribution")
acf(A_burn,main = "Autocorrelation function of join probability distribution")
#diagnostic for P
#-----


acceptance_rate_p  = (acc.count_p/(N_iter*K_true))*100
print(acceptance_rate_p)

plot(ts(alpha_burn), main = "Traceplot of alpha values")
abline(h = alpha, col = "red", lty = 2)

acf(alpha_burn,main = "Autocorrelation plot of alpha values")

mean_alpha= mean(alpha_burn)
print(mean_alpha)

alpha_df= data.frame(y=alpha_burn)

sd_alpha = sqrt(var(alpha_burn))

ggplot(alpha_df, aes(y)) +
  geom_density(fill = "purple", alpha = 0.5) +
  geom_vline(aes(xintercept = alpha, linetype = "True Value"), col = "red") +
  geom_vline(aes(xintercept = mean_alpha, linetype = "Mean Value"), col = "blue") +
  labs(x = "Density of alpha", y = "alpha values")

ggplot(alpha_df, aes(y)) +
  geom_density(fill = "purple", alpha = 0.5) +
  geom_vline(aes(xintercept = alpha, linetype = "True Value"), col = "red") +
  geom_pointrange(aes(x = mean_alpha, y = 0), 
                  xmin = max(mean_alpha - 1.96 * sd_alpha,0),
                  xmax = mean_alpha + 1.96 * sd_alpha,
                  color = "black") +
  labs(x = "Density of alpha", y = "alpha values")

table = matrix(0,K,K)
for(i in 1:K_true){
  for(j in 1:K_true){
    table[i,j]=((mean(p.container[i,j,]) - synth$P_matrix[i,j])*100)
  }
}
table%>% pander::pander()

diagnostic_df = data.frame(t(diagnostic))
names(diagnostic_df) = c('y.n',"my_alpha", "err1","err2","exp(r)")


diagnostic_df=diagnostic_df[-N_iter,] 

ggplot(alpha_df, aes(y)) +
   geom_density(fill = "purple", alpha = 0.5) +
    geom_vline(aes(xintercept = alpha, linetype = "True Value"), col = "red") +
    geom_pointrange(aes(x = mean_alpha, y = 0), 
                                       xmin = max(mean_alpha - 1.96 * sd_alpha,0),
                                        xmax = mean_alpha + 1.96 * sd_alpha,
                                        color = "black") +
    labs(x = "Density of alpha", y = "alpha values")



ggplot(diagnostic_df, aes(x = my_alpha,y=err1, col=factor(y.n)))+
  geom_point()+
  geom_vline(aes(xintercept= alpha), col="red")

ggplot(diagnostic_df, aes(x = my_alpha,y=err2, col=factor(y.n)))+
  geom_point()+
  geom_vline(aes(xintercept= alpha), col="red")

ggplot(diagnostic_df, aes(x = my_alpha, col=factor(y.n)))+
  geom_histogram(aes(fill=factor(y.n)))+
  geom_vline(aes(xintercept= alpha), col="red")


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



















