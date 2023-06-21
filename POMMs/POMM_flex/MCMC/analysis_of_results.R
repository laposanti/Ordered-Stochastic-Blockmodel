
samples_simple$z_container
# Stop parallel computing
stopCluster(cl)

seed79141<-readRDS("ModelPOMMPOMMK4_overlap0.2_alpha0.5_seed79141.RDS")
ts.plot(seed79141$A_container[-c(1:10000)])
acf(seed79141$A_container[-c(1:10000)])

seed79141$p_container

matrixP_list <- list()
p_true_list <- list()



file_names <- list.files(pattern = "ModelPOMMPOMMK4_overlap0.2_alpha0.5_seed")

for (file_name in file_names) {
  obj <- readRDS(file_name)
  matrixP <- obj$p_container
  true_P <-obj$p_true
  matrixP_list[[file_name]] <- matrixP
  p_true_list[[file_name]] <- true_P
  print(true_P)
}

blue_shades <- generate_color_gradient_K(5)

#Diagnostics for the P entries ---------
entry1=1
entry2=4
#checking convergence
density_plot_df<-data.frame(chain1_density = matrixP_list$ModelPOMMPOMMK4_overlap0.2_alpha0.5_seed79141.RDS[entry1,entry2,],
                            chain2_density = matrixP_list$ModelPOMMPOMMK4_overlap0.2_alpha0.5_seed79142.RDS[entry1,entry2,],
                            chain3_density = matrixP_list$ModelPOMMPOMMK4_overlap0.2_alpha0.5_seed79143.RDS[entry1,entry2,],
                            chain4_density = matrixP_list$ModelPOMMPOMMK4_overlap0.2_alpha0.5_seed79144.RDS[entry1,entry2,],
                            chain4_density = matrixP_list$ModelPOMMPOMMK4_overlap0.2_alpha0.5_seed79145.RDS[entry1,entry2,])




# Create a density plot of points for each level
ggplot() +
  # Add a layer for each level
  lapply(seq(1,5), function(i) {
    geom_density(data = data.frame(x = density_plot_df[,i]), aes(x = x, y = ..density.., color = paste("Chain ", i)), alpha = .5)
  }) +
  geom_vline(aes(xintercept = p_true_list$ModelPOMMPOMMK4_overlap0.2_alpha0.5_seed79141.RDS[entry1,entry2]), col="red", alpha=0.5)+
  # Set the x-axis limits
  scale_x_continuous(limits = c(min(density_plot_df), max(density_plot_df))) +# Remove the legend title
  guides(color = guide_legend(title = NULL)) +
  # Set the legend title
  labs(fill = "Chain", x = "Points", y = "Density", title = paste("Density Plot of Entries","P[",entry1,",",entry2,"]"))+
  # Set blue color scale for the lines
  scale_color_manual(values=blue_shades)+
  theme_bw()




for (file_name in file_names) {
  obj <- readRDS(file_name)
  z_container <- samples_simple$z_container
  z_true <- obj$z_true
  similarity_matrix = pr_cc(z_container[,-c(1:20000)])
  point_est = minVI(similarity_matrix)$cl
  print(adj.rand.index(point_est, samples_simple$z_true))
  mean(overlap_container)
  A_container <- obj$A_container
  z_MAP= z_container[,which(A_container == max(A_container))]
  print(adj.rand.index(z_MAP, z_true))
  ts.plot(A_container[-c(1:N_iter*0.25)])
}




for (file_name in file_names) {
  obj <- readRDS(file_name)
  z0 <- obj$init$z0
  overlap0 <- obj$init$overlap0
  alpha0 <- obj$init$alpha0
  print(adj.rand.index(obj$z_true, z0))
  print(abs(overlap0- obj$overlap))
  print(abs(overlap0- obj$alpha))
}

ts.plot(samples_simple$A_container[-c(1:25000)])
acf(samples_simple$A_container[-c(1:25000)])










