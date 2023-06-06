
sample_beta_trunc = function(N,alpha,beta,a,b){
  u = runif(N)
  return( qbeta(((pbeta(b, alpha, beta) -pbeta(a, alpha, beta))*u+pbeta(a, alpha, beta)),alpha,beta))}


for (i in seq_along(c_values)) {
  c <- c_values[i]
  beta_0_sum <- matrix(NA, nrow=K, ncol=N)
  for (j in 1:N) {
    b <- sample_beta_trunc(K,1,5,0,0.25)
    p <- numeric(n)
    p[1] <- 0.5
    for (k in 2:K) {
      p[k] <- b[k] * prod((0.75-0.5)- b[1:(k-1)])
    }
    p_sum <-  cumsum(p)
    beta_0_sum[,j] <-  (beta_max - 0.5) * (p_sum - min(p_sum)) / (max(p_sum) - min(p_sum)) + 0.5
  }
  beta_0_avg[[i]] <- rowSums(beta_0_sum)/ N
}
########################################
#Trying different concentration values
########################################
K <- 10
c_values <- c(0.3,0.5, 1,2,3, 10)  # modify as needed
N <- 1000  # number of samples per c value
num_c <- length(c_values)
beta_0_avg <- vector("list", num_c)

for (i in seq_along(c_values)) {
  c <- c_values[i]
  c=2
  beta_0_sum <- matrix(NA, nrow=K, ncol=N)
  for (j in 1:N) {
    b <- rbeta(K, 1, c)
    p <- numeric(n)
    p[1] <- b[1]
    for (k in 2:K) {
      p[k] <- b[k] * prod(1 - b[1:(k-1)])
    }
    p_sum <-  cumsum(p)
    beta_0_sum[,j] <-  (beta_max - 0.5) * (p_sum - min(p_sum)) / (max(p_sum) - min(p_sum)) + 0.5
  }
  beta_0_avg[[i]] <- rowSums(beta_0_sum)/ N
}

# plot all beta_0 values on the same graph
colors <- rainbow(num_c)
plot(1:K, beta_0_avg[[1]], type = "l", ylim = range(unlist(beta_0_avg)), col =colors[1], main = "Simulation of beta with different c values",
     xlab = "Level Sets", ylab= "Probability", cex.main=.9, cex.lab = .9,cex.axis = .9 )
for (i in 2:num_c) {
  lines(1:K, beta_0_avg[[i]], col = colors[i])
}
legend("bottomright", legend = paste0("c = ", c_values), lty = 1, col = colors, cex = 0.5)
#################################
#Trying different beta_max values
#################################
K <- 10
c= 3
beta_max_values <- c(0.6, 0.7, 0.8, 0.9)  # modify as needed
N <- 1000  # number of samples per c value
num_beta_max <- length(beta_max_values)
beta_0_avg <- vector("list", num_beta_max)

for (i in seq_along(beta_max_values)) {
  beta_max <- beta_max_values[i]
  beta_0_sum <- matrix(NA, nrow=K, ncol=N)
  for (j in 1:N) {
    b <- rbeta(K, 1, c)
    p <- numeric(K)
    p[1] <- b[1] 
    for (k in 2:K) {
      p[k] <- b[k] * prod(1 - b[1:(k-1)])
    }
    p_sum <-  cumsum(p)
    beta_0_sum[,j] <-  (beta_max - 0.5) * (p_sum - min(p_sum)) / (max(p_sum) - min(p_sum)) + 0.5
  }
  beta_0_avg[[i]] <- rowSums(beta_0_sum)/ N
}

# plot all beta_0 values on the same graph
colors <- rainbow(num_beta_max)
plot(1:K, beta_0_avg[[1]], type = "l", ylim = range(unlist(beta_0_avg)), col =colors[1], main = "Simulation of beta with different beta_max values",
     xlab = "Level Sets", ylab= "Probability", cex.main=.9, cex.lab = .9,cex.axis = .9 )
for (i in 2:num_beta_max) {
  lines(1:K, beta_0_avg[[i]], col = colors[i])
}
legend("topleft", legend = paste0("beta_max = ", beta_max_values), lty = 1, col = colors, cex = 0.5)


################################################
#Generating a matrix according to this new model
###############################################

simulating_POMM_stick_breaking = function(K,  c = .5, beta_max){
  
  #simulating from the prior via a stick-breaking process
  #simulating K beta values with a concentration parameter c
  b <- rbeta(K, 1, c)
  p <- numeric(K)
  p[1] <- 0.5
  # stick breaking prior on a stick of length 1
  for(i in 2:K){
    p[i] <-  b[i] * prod(1 - b[1:(i-1)])}
  #we sum over the sampled proportions so to ensure an increasing behaviour
  p_sum <-  cumsum(p)
  #we map the proportions on a desired scale
  beta_0 <-  (beta_max - 0.5) * (p_sum - min(p_sum)) / (max(p_sum) - min(p_sum)) + 0.5
  #L_k is a list containing all the diagonals of the upper triangular matrix
  L_k = diag_split_NA(K)
  #the diagonal is made by 0.5 values
  L_k[[K]] = rep(0.5,K)
  # Compute beta values
   for(i in (K-1):1){
    n_items = length(L_k[[i]])
    #generating the beta values using the proportions as truncations
    L_k[[i]] = sample_beta_trunc(n_items,1,1,beta_0[K-i], beta_0[K-i+1])
  }

  #putting back everything together in the matrix
  result = diag_matrix(K)
  for(i in 0:(K-1)){
    result[which(result==-i)]<-L_k[[K-i]]
  }
  result = semi_symmetric(as.matrix(result))
  return(result)}

###############################################
#Plotting the heatmap of a single configuration
###############################################
test = simulating_POMM_stick_breaking(10,c=3,beta_max = 0.8)

# convert matrix to data frame
df <- as.data.frame(t(test))
x = c(1:K)
y= c(K:1)
data <- expand.grid(X=x, Y=y)
df_long <- reshape2::melt(df)

df_plot = data.frame(X =data$X, Y= data$Y, Z = df_long$value)
# Create heatmap
ggplot(df_plot, aes(x=X, y=Y, fill=Z)) +
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  scale_x_discrete(limits=as.character(x), name="Column")+
  scale_y_discrete(limits=(as.character(y)), name="Row")+
  labs(title="Heatmap of POMM-generated matrix", x="Column", y="Row", fill="Value")+
  geom_text(aes(label = round(Z, 2)), size =2)




