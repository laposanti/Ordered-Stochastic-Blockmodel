
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

#alpha must be greater than zero.
#alpha in [0,1] --> truncations are a concave increasing function
#alpha ==1 --> truncations are a linearly increasing function
#alpha >1--> truncations are a convex increasing function
levels_alpha = list("concave" = .5, "linear"=1, "convex"=1.5)
#alpha = sample_norm_trunc(N = 1,m = levels_alpha$concave, s = var_alpha, a = 0,b=Inf)
#beta_max is defined on [0.5,1]. It controls the highest possible value in the matrix

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




grad_colors <- colorRampPalette(c("blue", "purple"))

# Density plot for the level sets with gradient colors
# Density plot for the level sets with gradient colors
density_concave_plot <- ggplot() +
  # Add a layer for each level with gradient fill colors
  lapply(seq_along(level_list_concave), function(i) {
    geom_density(data = data.frame(x = level_list_concave[[i]], i = i), 
                 aes(x = x, y = ..density.., fill = i), 
                 alpha = 0.5)
  }) +
  # Set the x-axis limits
  scale_x_continuous(limits = c(0.5, beta_max)) +
  # Set the legend title and gradient colors
  labs(fill = "Levels", x = "Points", y = "Density", title = "Concave case: Density Plot for Each Level") +
  scale_fill_gradientn(colors = grad_colors(5))


#heatmap of the generate SST matrix plot #################
heatmap_concave_plot = heat_map_blue(mean_concave_m,"Concave case: SST matrix using Power Law Prior" )

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


# create data frame for plotting
concave_df <- data.frame(x = 1:K,
                        y = empirical_concave_t,
                        mid = (concave_quantile95+concave_quantile5)/2,
                        y5 = concave_quantile5,
                        y95 =concave_quantile95 )

# Confidence interval plot ###################
confidence_concave_plot = ggplot(concave_df, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y5, ymax = y95), alpha = 0.2, fill = "blue") +
  geom_line(aes(y = y)) +
  geom_line(aes(y = mid), linetype = "dashed") +
  scale_x_continuous(limits = c(1, nrow(truncation_df))) +
  xlab("Index") +
  ylab("Observation Value") +
  ggtitle("Observations within Truncations Interval") +
  theme_bw()


grid.arrange(heatmap_concave_plot, confidence_concave_plot, density_concave_plot, ncol = 2, nrow=2,top = paste("Var(alpha) =", var_alpha,"|| beta_max =", beta_max,"|| K =", K ))


#########
#Linear case
######

#linear
for(j in 1:N){
  alpha = sample_norm_trunc(N = 1,m = levels_alpha$linear, s = var_alpha, a = .90,b=1.10)
  linear_matrix[,,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$matrix
  linear_truncations[,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$truncations
}

# Combine the four levels into a list
level_list_linear <- generalized_levels(linear_matrix,K=5,N = N)

grad_colors <- colorRampPalette(c("blue", "purple"))

# Density plot for the level sets with gradient colors
# Density plot for the level sets with gradient colors
density_linear_plot <- ggplot() +
  # Add a layer for each level with gradient fill colors
  lapply(seq_along(level_list_linear), function(i) {
    geom_density(data = data.frame(x = level_list_linear[[i]], i = i), 
                 aes(x = x, y = ..density.., fill = i), 
                 alpha = 0.5)
  }) +
  # Set the x-axis limits
  scale_x_continuous(limits = c(0.5, beta_max)) +
  # Set the legend title and gradient colors
  labs(fill = "Levels", x = "Points", y = "Density", title = "Linear case: Density Plot for Each Level") +
  scale_fill_gradientn(colors = grad_colors(5))

#linear mean
mean_linear_m = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    mean_linear_m[i,j] = sum(linear_matrix[i,j,])/N
  }}

mean_linear_t = rowSums(linear_truncations)/N
#heatmap of the generate SST matrix plot #################
heatmap_linear_plot = heat_map_blue(mean_linear_m,"Linear case: SST matrix using Power Law Prior" )


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

# create data frame for plotting
linear_df <- data.frame(x = 1:K,
                 y = empirical_linear_t,
                 mid = (linear_quantile95+linear_quantile5)/2,
                 y5 = linear_quantile5,
                 y95 =linear_quantile95 )

confidence_linear_plot = ggplot(linear_df, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y5, ymax = y95), alpha = 0.2, fill = "blue") +
  geom_line(aes(y = y)) +
  geom_line(aes(y = mid), linetype = "dashed") +
  scale_x_continuous(limits = c(1, K)) +
  xlab("Index") +
  ylab("Observation Value") +
  ggtitle("Observations within Truncations Interval") +
  theme_bw()

grid.arrange(heatmap_linear_plot, confidence_linear_plot, density_linear_plot, ncol = 2, nrow=2,top = paste("Var(alpha) =", var_alpha,"|| beta_max =", beta_max,"|| K =", K ))







##############
#Convex mean
#############


for(j in 1:N){
  alpha= alpha = sample_norm_trunc(N = 1,m = levels_alpha$convex, s = var_alpha, a = 1,b=Inf)
  convex_matrix[,,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$matrix
  convex_truncations[,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$truncations
}

# Combine the four levels into a list
level_list_convex <- generalized_levels(convex_matrix,K=5,N = N)


grad_colors <- colorRampPalette(c("blue", "purple"))

# Density plot for the level sets with gradient colors
# Density plot for the level sets with gradient colors
density_convex_plot <- ggplot() +
  # Add a layer for each level with gradient fill colors
  lapply(seq_along(level_list_convex), function(i) {
    geom_density(data = data.frame(x = level_list_convex[[i]], i = i), 
                 aes(x = x, y = ..density.., fill = i), 
                 alpha = 0.5)
  }) +
  # Set the x-axis limits
  scale_x_continuous(limits = c(0.5, beta_max)) +
  # Set the legend title and gradient colors
  labs(fill = "Levels", x = "Points", y = "Density", title = "Convex case: Density Plot for Each Level") +
  scale_fill_gradientn(colors = grad_colors(5))

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

#heatmap of the generate SST matrix plot #################
heatmap_convex_plot = heat_map_blue(mean_convex_m,"Convex case: SST matrix using Power Law Prior" )


empirical_convex_t = empirical_trunc(mean_convex_m, fun = "mean")
quantile5_convex = empirical_trunc(quantile_convex_m5, fun = "min")
quantile95_convex = empirical_trunc(quantile_convex_m95, fun = "max")

# truncation_df = matrix(0,K,3)
# truncation_df[1,] =.5
# for(i in 2:K){
#   truncation_df[i,1] = empirical_convex_t[i]
#   truncation_df[i,2]= mean_convex_t[i-1]
#   truncation_df[i,3]= mean_convex_t[i]
# }

# create data frame for plotting
convex_df <- data.frame(x = c(1:K),
                 y = empirical_convex_t,
                 mid = (quantile95_convex+quantile5_convex)/2,
                 y5 = quantile5_convex,
                 y95 =quantile95_convex )

# create plot
confidence_convex_plot = ggplot(convex_df, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y5, ymax = y95), alpha = 0.2, fill = "blue") +
  geom_line(aes(y = y)) +
  geom_line(aes(y = mid), linetype = "dashed") +
  scale_x_continuous(limits = c(1, K)) +
  xlab("Index") +
  ylab("Observation Value") +
  ggtitle("Observations within Truncations Interval") +
  theme_bw()



grid.arrange(heatmap_convex_plot, confidence_convex_plot, density_convex_plot, ncol = 2, nrow=2,top = paste("Var(alpha) =", var_alpha,"|| beta_max =", beta_max,"|| K =", K ))

#################
# Studying mean behaviour
#####


# Simulation Study
set.seed(123)

#setting the number of rows of the square matrix

#alpha must be greater than zero.
#alpha in [0,1] --> truncations are a concave increasing function
#alpha ==1 --> truncations are a linearly increasing function
#alpha >1--> truncations are a convex increasing function
levels_alpha_mean = list("medium_low" = .02, "low" = .05, "medium_low"=.1, "medium"=0.2)
#alpha = sample_norm_trunc(N = 1,m = levels_alpha$concave, s = var_alpha, a = 0,b=Inf)
#beta_max is defined on [0.5,1]. It controls the highest possible value in the matrix

N= 10000
K=5
var_alpha=.1
beta_max = .75

extremely_low_matrix = array(0, dim = c(K,K,N))
extremely_low_truncations = matrix(0, (K),N)

####
#extremely_low mean
####

for(j in 1:N){
  alpha = sample_norm_trunc(N = 1,m = levels_alpha_mean$extremely_low, s = var_alpha, a = 0.01,b=1)
  extremely_low_matrix[,,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$matrix
  extremely_low_truncations[,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$truncations
}


# Combine the four levels into a list
level_list_extremely_low <- generalized_levels(extremely_low_matrix,K=5,N = N)



density_extremely_low_plot <- ggplot() +
  # Add a layer for each level with gradient fill colors
  lapply(seq_along(level_list_extremely_low), function(i) {
    geom_density(data = data.frame(x = level_list_extremely_low[[i]], i = i), 
                 aes(x = x, y = ..density.., fill = i), 
                 alpha = 0.5)
  }) +
  # Set the x-axis limits
  scale_x_continuous(limits = c(0.5, beta_max)) +
  # Set the legend title and gradient colors
  labs(fill = "Levels", x = "Points", y = "Density", title = "Convex case: Density Plot for Each Level") +
  scale_fill_gradientn(colors = grad_colors(5))

mean_extremely_low_matrix_m = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    mean_extremely_low_matrix_m[i,j] = mean(extremely_low_matrix[i,j,])
    
  }}

#convex 5% quantile
quantile_extremely_low_matrix_m5 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_extremely_low_matrix_m5[i,j] =quantile(extremely_low_matrix[i,j,],probs = .05)
  }}
#convex 95% quantile
quantile_extremely_low_matrix_m95 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_extremely_low_matrix_m95[i,j] =quantile(extremely_low_matrix[i,j,],probs = .95)
  }}

#heatmap of the generate SST matrix plot #################
heatmap_convex_plot = heat_map_blue(mean_extremely_low_matrix_m,"Convex case: SST matrix using Power Law Prior" )


empirical_extremely_low_t = empirical_trunc(mean_extremely_low_matrix_m, fun = "mean")
quantile5_extremely_low = empirical_trunc(quantile_extremely_low_matrix_m5, fun = "min")
quantile95_extremely_low = empirical_trunc(quantile_extremely_low_matrix_m95, fun = "max")

# truncation_df = matrix(0,K,3)
# truncation_df[1,] =.5
# for(i in 2:K){
#   truncation_df[i,1] = empirical_convex_t[i]
#   truncation_df[i,2]= mean_convex_t[i-1]
#   truncation_df[i,3]= mean_convex_t[i]
# }

# create data frame for plotting
extremely_low_df <- data.frame(x = c(1:K),
                        y = empirical_extremely_low_t,
                        mid = (quantile95_extremely_low+quantile5_extremely_low)/2,
                        y5 = quantile5_extremely_low,
                        y95 =quantile95_extremely_low )

# create plot
confidence_convex_plot = ggplot(extremely_low_df, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y5, ymax = y95), alpha = 0.2, fill = "blue") +
  geom_line(aes(y = y)) +
  geom_line(aes(y = mid), linetype = "dashed") +
  scale_x_continuous(limits = c(1, K)) +
  xlab("Index") +
  ylab("Observation Value") +
  ggtitle("Observations within Truncations Interval") +
  theme_bw()



grid.arrange(heatmap_convex_plot, confidence_convex_plot, density_convex_plot, ncol = 2, nrow=2,top = paste("Var(alpha) =", var_alpha,"|| beta_max =", beta_max,"|| K =", K ))


####
#low mean
####

for(j in 1:N){
  alpha = sample_norm_trunc(N = 1,m = levels_alpha_mean$low, s = var_alpha, a = 0.01,b=1)
  low_matrix[,,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$matrix
  low_truncations[,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$truncations
}


# Combine the four levels into a list
level_list_low <- generalized_levels(low_matrix,K=5,N = N)



density_low_plot <- ggplot() +
  # Add a layer for each level with gradient fill colors
  lapply(seq_along(level_list_low), function(i) {
    geom_density(data = data.frame(x = level_list_low[[i]], i = i), 
                 aes(x = x, y = ..density.., fill = i), 
                 alpha = 0.5)
  }) +
  # Set the x-axis limits
  scale_x_continuous(limits = c(0.5, beta_max)) +
  # Set the legend title and gradient colors
  labs(fill = "Levels", x = "Points", y = "Density", title = "Convex case: Density Plot for Each Level") +
  scale_fill_gradientn(colors = grad_colors(5))

mean_low_matrix_m = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    mean_low_matrix_m[i,j] = mean(low_matrix[i,j,])
    
  }}

#convex 5% quantile
quantile_low_matrix_m5 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_low_matrix_m5[i,j] =quantile(low_matrix[i,j,],probs = .05)
  }}
#convex 95% quantile
quantile_low_matrix_m95 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_low_matrix_m95[i,j] =quantile(low_matrix[i,j,],probs = .95)
  }}

#heatmap of the generate SST matrix plot #################
heatmap_convex_plot = heat_map_blue(mean_low_matrix_m,"Convex case: SST matrix using Power Law Prior" )


empirical_low_t = empirical_trunc(mean_low_matrix_m, fun = "mean")
quantile5_low = empirical_trunc(quantile_low_matrix_m5, fun = "min")
quantile95_low = empirical_trunc(quantile_low_matrix_m95, fun = "max")

# truncation_df = matrix(0,K,3)
# truncation_df[1,] =.5
# for(i in 2:K){
#   truncation_df[i,1] = empirical_convex_t[i]
#   truncation_df[i,2]= mean_convex_t[i-1]
#   truncation_df[i,3]= mean_convex_t[i]
# }

# create data frame for plotting
low_df <- data.frame(x = c(1:K),
                               y = empirical_low_t,
                               mid = (quantile95_low+quantile5_low)/2,
                               y5 = quantile5_low,
                               y95 =quantile95_low )

# create plot
confidence_convex_plot = ggplot(low_df, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y5, ymax = y95), alpha = 0.2, fill = "blue") +
  geom_line(aes(y = y)) +
  geom_line(aes(y = mid), linetype = "dashed") +
  scale_x_continuous(limits = c(1, K)) +
  xlab("Index") +
  ylab("Observation Value") +
  ggtitle("Observations within Truncations Interval") +
  theme_bw()



grid.arrange(heatmap_convex_plot, confidence_convex_plot, density_convex_plot, ncol = 2, nrow=2,top = paste("Var(alpha) =", var_alpha,"|| beta_max =", beta_max,"|| K =", K ))

####
#medium mean
####

for(j in 1:N){
  alpha = sample_norm_trunc(N = 1,m = levels_alpha_mean$medium, s = var_alpha, a = 0.01,b=1)
  medium_matrix[,,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$matrix
  medium_truncations[,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$truncations
}


# Combine the four levels into a list
level_list_medium <- generalized_levels(medium_matrix,K=5,N = N)



density_medium_plot <- ggplot() +
  # Add a layer for each level with gradient fill colors
  lapply(seq_along(level_list_medium), function(i) {
    geom_density(data = data.frame(x = level_list_medium[[i]], i = i), 
                 aes(x = x, y = ..density.., fill = i), 
                 alpha = 0.5)
  }) +
  # Set the x-axis limits
  scale_x_continuous(limits = c(0.5, beta_max)) +
  # Set the legend title and gradient colors
  labs(fill = "Levels", x = "Points", y = "Density", title = "Convex case: Density Plot for Each Level") +
  scale_fill_gradientn(colors = grad_colors(5))

mean_medium_matrix_m = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    mean_medium_matrix_m[i,j] = mean(medium_matrix[i,j,])
    
  }}

#convex 5% quantile
quantile_medium_matrix_m5 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_medium_matrix_m5[i,j] =quantile(medium_matrix[i,j,],probs = .05)
  }}
#convex 95% quantile
quantile_medium_matrix_m95 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_medium_matrix_m95[i,j] =quantile(medium_matrix[i,j,],probs = .95)
  }}

#heatmap of the generate SST matrix plot #################
heatmap_convex_plot = heat_map_blue(mean_medium_matrix_m,"Convex case: SST matrix using Power Law Prior" )


empirical_medium_t = empirical_trunc(mean_medium_matrix_m, fun = "mean")
quantile5_medium = empirical_trunc(quantile_medium_matrix_m5, fun = "min")
quantile95_medium = empirical_trunc(quantile_medium_matrix_m95, fun = "max")

# truncation_df = matrix(0,K,3)
# truncation_df[1,] =.5
# for(i in 2:K){
#   truncation_df[i,1] = empirical_convex_t[i]
#   truncation_df[i,2]= mean_convex_t[i-1]
#   truncation_df[i,3]= mean_convex_t[i]
# }

# create data frame for plotting
medium_df <- data.frame(x = c(1:K),
                               y = empirical_medium_t,
                               mid = (quantile95_medium+quantile5_medium)/2,
                               y5 = quantile5_medium,
                               y95 =quantile95_medium )

# create plot
confidence_convex_plot = ggplot(medium_df, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y5, ymax = y95), alpha = 0.2, fill = "blue") +
  geom_line(aes(y = y)) +
  geom_line(aes(y = mid), linetype = "dashed") +
  scale_x_continuous(limits = c(1, K)) +
  xlab("Index") +
  ylab("Observation Value") +
  ggtitle("Observations within Truncations Interval") +
  theme_bw()



grid.arrange(heatmap_convex_plot, confidence_convex_plot, density_convex_plot, ncol = 2, nrow=2,top = paste("Var(alpha) =", var_alpha,"|| beta_max =", beta_max,"|| K =", K ))


####
#medium_low mean
####

for(j in 1:N){
  alpha = sample_norm_trunc(N = 1,m = levels_alpha_mean$medium_low, s = var_alpha, a = 0.01,b=1)
  medium_low_matrix[,,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$matrix
  medium_low_truncations[,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$truncations
}


# Combine the four levels into a list
level_list_medium_low <- generalized_levels(medium_low_matrix,K=5,N = N)



density_medium_low_plot <- ggplot() +
  # Add a layer for each level with gradient fill colors
  lapply(seq_along(level_list_medium_low), function(i) {
    geom_density(data = data.frame(x = level_list_medium_low[[i]], i = i), 
                 aes(x = x, y = ..density.., fill = i), 
                 alpha = 0.5)
  }) +
  # Set the x-axis limits
  scale_x_continuous(limits = c(0.5, beta_max)) +
  # Set the legend title and gradient colors
  labs(fill = "Levels", x = "Points", y = "Density", title = "Convex case: Density Plot for Each Level") +
  scale_fill_gradientn(colors = grad_colors(5))

mean_medium_low_matrix_m = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    mean_medium_low_matrix_m[i,j] = mean(medium_low_matrix[i,j,])
    
  }}

#convex 5% quantile
quantile_medium_low_matrix_m5 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_medium_low_matrix_m5[i,j] =quantile(medium_low_matrix[i,j,],probs = .05)
  }}
#convex 95% quantile
quantile_medium_low_matrix_m95 = matrix(0,K,K)
for(i in 1:K){
  for(j in 1:K){
    quantile_medium_low_matrix_m95[i,j] =quantile(medium_low_matrix[i,j,],probs = .95)
  }}

#heatmap of the generate SST matrix plot #################
heatmap_convex_plot = heat_map_blue(mean_medium_low_matrix_m,"Convex case: SST matrix using Power Law Prior" )


empirical_medium_low_t = empirical_trunc(mean_medium_low_matrix_m, fun = "mean")
quantile5_medium_low = empirical_trunc(quantile_medium_low_matrix_m5, fun = "min")
quantile95_medium_low = empirical_trunc(quantile_medium_low_matrix_m95, fun = "max")

# truncation_df = matrix(0,K,3)
# truncation_df[1,] =.5
# for(i in 2:K){
#   truncation_df[i,1] = empirical_convex_t[i]
#   truncation_df[i,2]= mean_convex_t[i-1]
#   truncation_df[i,3]= mean_convex_t[i]
# }

# create data frame for plotting
medium_low_df <- data.frame(x = c(1:K),
                               y = empirical_medium_low_t,
                               mid = (quantile95_medium_low+quantile5_medium_low)/2,
                               y5 = quantile5_medium_low,
                               y95 =quantile95_medium_low )

# create plot
confidence_convex_plot = ggplot(medium_low_df, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y5, ymax = y95), alpha = 0.2, fill = "blue") +
  geom_line(aes(y = y)) +
  geom_line(aes(y = mid), linetype = "dashed") +
  scale_x_continuous(limits = c(1, K)) +
  xlab("Index") +
  ylab("Observation Value") +
  ggtitle("Observations within Truncations Interval") +
  theme_bw()



grid.arrange(heatmap_convex_plot, confidence_convex_plot, density_convex_plot, ncol = 2, nrow=2,top = paste("Var(alpha) =", var_alpha,"|| beta_max =", beta_max,"|| K =", K ))


