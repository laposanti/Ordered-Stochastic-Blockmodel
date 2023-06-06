


N= 10000
K=5
var_alpha=2
beta_max = .75
mu = 4
extremely_low_matrix = array(0, dim = c(K,K,N))
extremely_low_truncations = matrix(0, (K),N)

####
#extremely_low mean
####

for(j in 1:N){
  alpha = sample_norm_trunc(N = 1,m = mu, s = var_alpha, a = 0.01,b=20)
  extremely_low_matrix[,,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$matrix
  extremely_low_truncations[,j] = simulating_POMM_powerlaw(K,alpha = alpha,beta_max = beta_max)$truncations
}


# Combine the four levels into a list
level_list_extremely_low <- generalized_levels(extremely_low_matrix,K=K,N = N)



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
  labs(fill = "Levels", x = "Points", y = "Density", title = "Density Plot for Each Level") +
  scale_fill_gradientn(colors = grad_colors(K-1))


mean_extremely_low_matrix_t =  matrix(0,K,1)
for(i in 1:K){
    mean_extremely_low_matrix_t[i] = mean(extremely_low_truncations[i,])
}

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
heatmap_convex_plot = heat_map_blue(mean_extremely_low_matrix_m,"SST matrix using Power Law Prior" )


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
                               true_u= mean_extremely_low_matrix_t,
                               true_l = c(0.5, mean_extremely_low_matrix_t[1:K-1]),
                               y5 = quantile5_extremely_low,
                               y95 =quantile95_extremely_low )

# create plot
confidence_convex_plot = ggplot(extremely_low_df, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = y5, ymax = y95), alpha = 0.2, fill = "blue") +
  geom_line(aes(y = y)) +
  geom_line(aes(y = mid), linetype = "dashed") +
  geom_line(aes(y = true_u), linetype = "dashed", col="red") +
  geom_line(aes(y = true_l), linetype = "dashed", col="red") +
  scale_x_continuous(limits = c(1, K)) +
  xlab("Index") +
  ylab("Observation Value") +
  ggtitle("Observations within Truncations Interval") +
  theme_bw()

plot_title = paste("mu_", mu,"_sigma_", var_alpha,"beta_max_", beta_max,"_K_", K )

plot_title= gsub(".", "_",as.character(plot_title),fixed = T)
plot_title= gsub(" ", "",as.character(plot_title))
plot_title= paste(plot_title,".png")
plot_title= gsub(" ", "",as.character(plot_title))


setwd("/Users/lapo_santi/Desktop/Nial/project/POMMs/power-law prior/simulation_graphs/mu_study")
g = grid.arrange(heatmap_convex_plot, confidence_convex_plot, density_extremely_low_plot, ncol = 2, nrow=2,top = paste("mu= ", mu,"|| Var(alpha) =", var_alpha,"|| beta_max =", beta_max,"|| K =", K ))

ggsave(file=plot_title, g,width = 3000,height = 2000,units = "px")

rm(list = ls())

