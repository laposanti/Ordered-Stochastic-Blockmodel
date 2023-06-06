
empirical_trunc = function(matrix, fun){
  empirical_trunc = vector()
  for( i in K:1){
    if(fun == "mean"){
    empirical_trunc=append(empirical_trunc,mean(diag_split_matrix(matrix)[[i]]))
    }else if(fun == "min"){
      empirical_trunc=append(empirical_trunc,min(diag_split_matrix(matrix)[[i]]))
  }else if(fun =="max"){
    empirical_trunc=append(empirical_trunc,max(diag_split_matrix(matrix)[[i]]))
  }}
  return(empirical_trunc)}

generalized_levels <- function(concave_matrix, K, N) {
  K=K-1
  level_list <- vector("list", K)
  
  for (i in 1:N) {
    split_matrices <- diag_split_matrix(concave_matrix[,,i])
    
    for (j in 1:K) {
      level_list[[j]] <- append(level_list[[j]], split_matrices[[j]])
    }
  }
  
  return(level_list)
}



# Simulation Study
set.seed(123)

#setting the number of rows of the square matrix
K =10
#alpha must be greater than zero.
#alpha in [0,1] --> truncations are a concave increasing function
#alpha ==1 --> truncations are a linearly increasing function
#alpha >1--> truncations are a convex increasing function
levels_alpha = list("concave" = .5, "linear"=1, "convex"=1.4)
var_alpha  = 1
alpha = sample_norm_trunc(N = 1,m = levels_alpha$concave, s = var_alpha, a = 0,b=Inf)
#beta_max is defined on [0.5,1]. It controls the highest possible value in the matrix
beta_max = .75

ex= list(a = list(matrix(1,2,2),2),
     b=list(3,4),
     c=list(5,6))





N= 10000
K=5
var_alpha=.5
beta_max = .75
concave_matrix = array(0, dim = c(K,K,N))
concave_truncations = matrix(0, (K),N)
linear_matrix = array(0, dim = c(K,K,N))
linear_truncations = matrix(0, (K),N)
convex_matrix = array(0, dim = c(K,K,N))
convex_truncations = matrix(0, (K),N)

########
#Concave case
########

for(j in 1:N){
  alpha= alpha = sample_norm_trunc(N = 1,m = levels_alpha$concave, s = var_alpha, a = 0.01,b=1)
  concave_matrix[,,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$matrix
  concave_truncations[,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$truncations
}


# Combine the four levels into a list
level_list_concave <- generalized_levels(concave_matrix,K=5,N = N)


# Create a density plot of points for each level
ggplot() +
  # Add a layer for each level
  lapply(seq_along(level_list_concave), function(i) {
    geom_density(data = data.frame(x = level_list_concave[[i]]), aes(x = x, y = ..density.., fill = paste0("Level ", i)), alpha = 0.5)
  }) +
  # Set the x-axis limits
  scale_x_continuous(limits = c(0.5, beta_max)) +
  # Set the legend title
  labs(fill = "Levels", x = "Points", y = "Density", title = "Density Plot of Points for Each Level")



mean_concave_t = rowSums(concave_truncations)/N

plot(ts(mean_concave_t), col= 'blue')
lines(ts(empirical_trunc(mean_concave_m)), col = 'red')

heat_map_blue(mean_concave_m)

#concave mean
mean_concave_m = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    mean_concave_m[i,j] = mean(concave_matrix[i,j,])
  }}

#concave 5% quantile
quantile_concave_m5 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_concave_m5[i,j] =quantile(concave_matrix[i,j,],probs = .05)
  }}

#concave 95% quantile
quantile_concave_m95 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_concave_m95[i,j] =quantile(concave_matrix[i,j,],probs = .95)
  }}


empirical_concave_t = empirical_trunc(mean_concave_m,fun = "mean")
concave_quantile5 = empirical_trunc(quantile_concave_m5, fun = "min")
concave_quantile95 = empirical_trunc(quantile_concave_m95, fun ="max")

# truncation_df = matrix(0,K,3)
# truncation_df[1,] =.5
# for(i in 2:K){
#   truncation_df[i,1] = empirical_convex_t[i]
#   truncation_df[i,2]= mean_convex_t[i-1]
#   truncation_df[i,3]= mean_convex_t[i]
# }

# create data frame for plotting
concave_df <- data.frame(x = 1:nrow(truncation_df),
                        y = empirical_concave_t,
                        mid = (concave_quantile95+concave_quantile5)/2,
                        y5 = concave_quantile5,
                        y95 =concave_quantile95 )

# create plot
ggplot(concave_df, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y5, ymax = y95), alpha = 0.2, fill = "red") +
  geom_line(aes(y = y)) +
  geom_line(aes(y = mid), linetype = "dashed") +
  scale_x_continuous(limits = c(1, nrow(truncation_df))) +
  xlab("Index") +
  ylab("Observation Value") +
  ggtitle("Observations within Truncations Interval") +
  theme_bw()

#########
#Linear case
######

#linear
for(j in 1:N){
  alpha = sample_norm_trunc(N = 1,m = levels_alpha$concave, s = var_alpha, a = .90,b=1.10)
  linear_matrix[,,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$matrix
  linear_truncations[,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$truncations
}

# Combine the four levels into a list
level_list_linear <- generalized_levels(linear_matrix,K=5,N = N)


# Create a density plot of points for each level
ggplot() +
  # Add a layer for each level
  lapply(seq_along(level_list_linear), function(i) {
    geom_density(data = data.frame(x = level_list_linear[[i]]), aes(x = x, y = ..density.., fill = paste0("Level ", i)), alpha = 0.5)
  }) +
  # Set the x-axis limits
  scale_x_continuous(limits = c(0.5, beta_max)) +
  # Set the legend title
  labs(fill = "Levels", x = "Points", y = "Density", title = "Density Plot of Points for Each Level")

#linear mean
mean_linear_m = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    mean_linear_m[i,j] = sum(linear_matrix[i,j,])/N
}}
mean_linear_t = rowSums(linear_truncations)/N
heat_map_blue(mean_linear_m)
plot(ts(mean_linear_t), col= 'blue')
lines(ts(empirical_trunc(mean_linear_m)), col = 'red')

#linear 5% quantile
quantile_linear_m5 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_linear_m5[i,j] =quantile(linear_matrix[i,j,],probs = .05)
  }}
#linear 95% quantile
quantile_linear_m95 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_linear_m95[i,j] =quantile(linear_matrix[i,j,],probs = .95)
  }}


empirical_linear_t = empirical_trunc(mean_linear_m,fun = "mean")
linear_quantile5 = empirical_trunc(quantile_linear_m5,fun = "min")
linear_quantile95 = empirical_trunc(quantile_linear_m95,fun = "max")

# truncation_df = matrix(0,K,3)
# truncation_df[1,] =.5
# for(i in 2:K){
#   truncation_df[i,1] = empirical_convex_t[i]
#   truncation_df[i,2]= mean_convex_t[i-1]
#   truncation_df[i,3]= mean_convex_t[i]
# }

# create data frame for plotting
linear_df <- data.frame(x = 1:nrow(truncation_df),
                 y = empirical_linear_t,
                 mid = (linear_quantile95+linear_quantile5)/2,
                 y5 = linear_quantile5,
                 y95 =linear_quantile95 )

# create plot
ggplot(linear_df, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y5, ymax = y95), alpha = 0.2, fill = "red") +
  geom_line(aes(y = y)) +
  geom_line(aes(y = mid), linetype = "dashed") +
  scale_x_continuous(limits = c(1, nrow(truncation_df))) +
  xlab("Index") +
  ylab("Observation Value") +
  ggtitle("Observations within Truncations Interval") +
  theme_bw()




##############
#Convex mean
#############

#concave
for(j in 1:N){
  alpha= alpha = sample_norm_trunc(N = 1,m = levels_alpha$concave, s = var_alpha, a = 1,b=5)
  convex_matrix[,,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$matrix
  convex_truncations[,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$truncations
}

# Combine the four levels into a list
level_list_convex <- generalized_levels(convex_matrix,K=5,N = N)


# Create a density plot of points for each level
ggplot() +
  # Add a layer for each level
  lapply(seq_along(level_list_convex), function(i) {
    geom_density(data = data.frame(x = level_list_convex[[i]]), aes(x = x, y = ..density.., fill = paste0("Level ", i)), alpha = 0.5)
  }) +
  # Set the x-axis limits
  scale_x_continuous(limits = c(0.5, beta_max)) +
  # Set the legend title
  labs(fill = "Levels", x = "Points", y = "Density", title = "Density Plot of Points for Each Level")


mean_convex_m = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    mean_convex_m[i,j] = sum(convex_matrix[i,j,])/N
    
  }}

#convex 5% quantile
quantile_convex_m5 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_convex_m5[i,j] =quantile(convex_matrix[i,j,],probs = .05)
  }}
#convex 95% quantile
quantile_convex_m95 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_convex_m95[i,j] =quantile(convex_matrix[i,j,],probs = .95)
  }}


empirical_convex_t = empirical_trunc(mean_convex_m, fun = "mean")
quantile5 = empirical_trunc(quantile_convex_m5, fun = "min")
quantile95 = empirical_trunc(quantile_convex_m95, fun = "max")

# truncation_df = matrix(0,K,3)
# truncation_df[1,] =.5
# for(i in 2:K){
#   truncation_df[i,1] = empirical_convex_t[i]
#   truncation_df[i,2]= mean_convex_t[i-1]
#   truncation_df[i,3]= mean_convex_t[i]
# }

# create data frame for plotting
df <- data.frame(x = 1:nrow(truncation_df),
                 y = truncation_df[,1],
                 mid = (quantile95+quantile5)/2,
                 y5 = quantile5,
                 y95 =quantile95 )

# create plot
ggplot(df, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y5, ymax = y95), alpha = 0.2, fill = "red") +
  geom_line(aes(y = y)) +
  geom_line(aes(y = mid), linetype = "dashed") +
  scale_x_continuous(limits = c(1, nrow(truncation_df))) +
  xlab("Index") +
  ylab("Observation Value") +
  ggtitle("Observations within Truncations Interval") +
  theme_bw()




