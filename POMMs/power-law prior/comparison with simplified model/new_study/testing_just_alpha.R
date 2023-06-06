
library(fossil)
library(dplyr)
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")

#Hyperparameters
N=100
a=1
b=1
K=3

M=10000
max_clust=3
min_clust= 3


N_iter = 20000
alpha=1


gamma_vec = rep(1,min_clust)
beta_max=.8

set.seed(16)
synth = simulating_tournament_new(N,alpha,beta_max,K, M = M,gamma_vec = gamma_vec,n_ij_max = 6,model = "POMM",diag0.5 = T)

K_true=K
z.true = synth$z_true
labels_available=1:K

n_ij_matrix = synth$n_ij_true
y_ij_matrix= synth$y_ij_true



#selecting just those values such that there is inside a non-0 entry
upper.tri_n_ij = upper.tri(n_ij_matrix)
non_negative_n_ij = which(upper.tri_n_ij & n_ij_matrix > 0, arr.ind = T)
p_logical_mat = non_negative_n_ij 


#retrieving the non-0 values
n_ij = n_ij_matrix[p_logical_mat]
y_ij = y_ij_matrix[p_logical_mat]

p.true=synth$P_matrix
p_ij_true= synth$p_ij_true[p_logical_mat]
#p_ij_true= synth$p_ij_true[non_negative_n_ij]


cor(p_ij_true ,y_ij)
table(p_ij_true[non_negative_n_ij]== p_ij_true[p_logical_mat])


maty.n=p_logical_mat
# Convert matrix to long format data.frame
df <- reshape2::melt(maty.n)
names(df) <- c("var1", "var2", "y.n")

# Create heatmap
ggplot(df, aes(x = var2, y = desc(var1), fill = factor(y.n))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("white", "red")) +
  labs(x = "", y = "") +
  theme_bw()


df_1 = data.frame(cbind(non_negative_n_ij, rep(T,3400)))
names(df_1) <- c("var1", "var2", "y.n")
# Create heatmap
ggplot(df_1, aes(x = var2, y = desc(var1), fill = factor(y.n))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("red")) +
  labs(x = "", y = "") +
  theme_bw()

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

truncations_current = improper_prior5(K,beta_max,alpha)

#Setting up quantities needed within computations


p_logical_mat = (upper.tri(p_n_current,diag = T) & n_ij_matrix>0)
p_ij_current = p_n_current[p_logical_mat]

A_seq[1] = sum(dbinom(y_ij, n_ij, p_ij_current, log=T)) +  l_like_p_ij(p.true, truncations_current) +
  dlnorm_param(alpha_prime) + ddirichlet_multinomial(N,K_true,n_k = n_k_current,my_alpha = gamma_vec)


#p_ij_current =  p_n_current[non_negative_n_ij]

#containers for the counts of accepted proposals
acc.count_z = 0
acc.count_p = 0

#setting time tracker
pb=txtProgressBar(min=1,max=N_iter)
j=2

diagnostic = matrix(0, nrow=5,ncol=N_iter)

#READY TO BOMB!
while (j < N_iter + 1) {
  setTxtProgressBar(pb, j)
  
  alpha_current=1

  
  #Update of the P matrix
  #----
  sigma_prime = 1
  
  #proposing a new alpha
  alpha_prime <- sample_norm_trunc(1, alpha_current,s =sigma_prime,a = 0.01,b = 3)
  truncations_prime <- improper_prior5(K,beta_max,alpha = alpha_prime)


  
  r= (sum(dbinom(y_ij, n_ij, p_ij_true, log=T))+   l_like_p_ij(p.true, truncations_prime) + dlnorm_param(alpha_prime))  -
    (sum(dbinom(y_ij, n_ij, p_ij_true, log = T))+  l_like_p_ij(p.true, truncations_current) + dlnorm_param(alpha_current))
  
  alpha_r = min(1, exp(r))
  u = runif(1)
  if(r == -Inf){}else if(u<alpha_r ){
    diagnostic[1,j-1] = 1
    diagnostic[2,j-1] = alpha_prime
    diagnostic[5,j-1] = exp(r)
    #counting number of accepted proposals
    acc.count_p = acc.count_p+1
    #updating quantities
    alpha_current = alpha_prime
    truncations_current=truncations_prime
  }else{
    diagnostic[1,j-1] = 0
    diagnostic[2,j-1] = alpha_prime
    diagnostic[5,j-1] = exp(r)
  }
  #storing results for diagnostics
  alpha.container[j]= alpha_current
  
  A_seq[j] = sum(dbinom(y_ij, n_ij, p_ij_current, log=T)) +  l_like_p_ij(p.true,truncations = truncations_current) +
    dlnorm_param(alpha_current) + ddirichlet_multinomial(N,K_true,n_k = n_k_current,my_alpha = gamma_vec)
  j=j+1
}


ts.plot(A_seq[-c(1:10000)])
acf(A_seq[-c(1:10000)])
#diagnostic for P
#-----

acceptance_rate_p  = (acc.count_p/N_iter)*100
print(acceptance_rate_p)

plot(ts(alpha.container), main = "Traceplot of alpha values")
abline(h = alpha, col = "red", lty = 2)
acf(alpha.container[-c(1:10000)],main = "ACF of alpha values")

mean(alpha.container[-c(1:10000)])
alpha_df = data.frame(alpha_values = alpha.container)

ggplot(alpha_df, aes(x= alpha_values))+
  geom_density(fill="purple")+
  geom_vline(aes(xintercept= alpha), col="red")
  
plot(density(alpha.container))
abline(v= alpha, col = "red", lty = 2)+
  geom_vline(aes(xintercept= alpha), col="red")

for(i in 1:K_true){
  for(j in 1:K_true){
    print(mean(p.container[i,j,]) - synth$P_matrix[i,j])
  }
}

diagnostic_df = data.frame(t(diagnostic))
names(diagnostic_df) = c('y.n',"my_alpha", "err1","err2","exp(r)")
diagnostic_df=diagnostic_df[-N_iter,] 

ggplot(diagnostic_df, aes(x = my_alpha,y=exp(r), col=factor(y.n)))+
  geom_point()+
  scale_y_log10()+
  geom_vline(aes(xintercept= alpha), col="red")




ggplot(diagnostic_df, aes(x = my_alpha,y=err1, col=factor(y.n)))+
  geom_point()+
  geom_vline(aes(xintercept= alpha), col="red")

ggplot(diagnostic_df, aes(x = my_alpha,y=err2, col=factor(y.n)))+
  geom_point()+
  geom_vline(aes(xintercept= alpha), col="red")

est.mean = diagnostic_df%>% select(my_alpha,y.n)%>%filter(y.n==1)%>% na.omit()%>% apply(2,mean) %>%as.vector()
ggplot(diagnostic_df, aes(x = my_alpha, col=factor(y.n)))+
  geom_histogram(aes(fill=factor(y.n)))+
  geom_vline(aes(xintercept= alpha), col="red")+
  geom_vline(aes(xintercept= est.mean[1]), col="blue")



















