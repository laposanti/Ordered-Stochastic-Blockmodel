library(gt)
library(gridExtra)
library(ggplot2)
library(truncnorm)
K_list = c(3,5,9)
K=K_list[2]
S_list = c(0.01,0.5,3)
alpha_list = c(1, 1, 1)
n_samples=10000
S = .21
beta_max = .8
alpha=.5
diag0.5=T
true_alpha<-alpha



#creating a sample of P matrices
p_container = array(0, dim=c(K,K,n_samples))
for(i in 1:n_samples){
  trunc = improper_prior5(K,alpha = alpha,diag0.5 = diag0.5, beta_max = beta_max )
  p_container[,,i] = simulating_overlapping_POMM_powerlaw_norm(K, alpha=alpha, S  = S, beta_max = beta_max,truncations = trunc,diag0.5 = diag0.5)
}

level_list_p_container<- generalized_levels(p_container,K,n_samples, diag0.5 = diag0.5)
# Combine the four levels into a list

p_container_simple = array(0, dim=c(K,K,n_samples))
for(i in 1:n_samples){
  p_container_simple[,,i] = matrix(runif(K*K, 0.5, beta_max), K, K)
}

level_list_p_container_simple<- generalized_levels(p_container_simple,K,n_samples, diag0.5 = diag0.5)

p_ij_lst<- list()

for (i in 1:(K-1)){
  for(j in (i+1):K){
    p_ij_lst[[paste0("entry",i,j)]]<-p_container[i,j,]
  }
}

my_p<-list()

for(i in  1:length(p_ij_lst)){
  i_th_el<- data.frame(y=p_ij_lst[[i]])
  my_p[[names(p_ij_lst[i])]]<- ggplot(i_th_el)+
    geom_density(aes(x=y))+
    theme_bw()+
    labs(title = names(p_ij_lst[i]))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}

my_table<-list()

for(i in  1:length(p_ij_lst)){
  my_data_i <- p_ij_lst[[i]]
  my_mean_i<-mean(my_data_i)%>% round(2)
  IQR_i<-IQR(my_data_i)%>% round(2)
  sd_i<-sd(my_data_i)%>% round(2)
  my_name_i<-names(p_ij_lst[i])
  i_th_df<- data.frame(my_name_i= c(my_mean_i,paste0("[",my_mean_i - IQR_i," - ", my_mean_i + IQR_i,"]"),sd_i))
  colnames(i_th_df) <- my_name_i 
  rownames(i_th_df)<- c('avg','IQR','sd')
  my_table[[names(p_ij_lst[i])]] <- tableGrob(i_th_df)
  }

sqrt(0.08)

hlay <- rbind(c(NA,1,2,3,4),
              c(11,NA,5,6,7),
              c(12,15,NA,8,9),
              c(13,16,18,NA,10),
              c(14,17,19,20,NA))

grid.arrange(
  my_p[[1]],
  my_p[[2]],
  my_p[[3]],
  my_p[[4]],
  my_p[[5]],
  my_p[[6]],
  my_p[[7]],my_p[[8]],my_p[[9]],my_p[[10]],
  my_table[[1]],
  my_table[[2]],
  my_table[[3]],
  my_table[[4]],
  my_table[[5]],
  my_table[[6]],
  my_table[[7]],my_table[[8]],my_table[[9]],my_table[[10]],
  layout_matrix = hlay,
  top=paste0("P_ij density for K= ",K," alpha= ",alpha," sigma= ",S))
 
















