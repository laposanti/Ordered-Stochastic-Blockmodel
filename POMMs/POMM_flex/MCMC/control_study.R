

source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/functions_container_flex.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
source("/Users/lapo_santi/Desktop/Nial/project/simplified model/SaraWade.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/adaptive_POMM_MCMC_function.R")
source("/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/POMM_flex/MCMC/Inference_functions.R")

library(EnvStats)
library(ggplot2)
library(dplyr)
library(ggside)
library(truncnorm)
library(fossil)
M= 13  
N = 30
N_iter =20000
N_ij = matrix(M,N,N)
alpha=.5
S=.01
K=13
beta_max=0.8
gamma_vec = rep(1/K,K)
diag0.5=T
trunc = improper_prior5(K,beta_max ,alpha,diag0.5)
P=simulating_overlapping_POMM_powerlaw_norm(K,alpha,S,trunc,beta_max,diag0.5)

z = rep_len(1:K, N)

z_mat=vec2mat(z)
P_NbyN<-calculate_victory_probabilities(z_mat,P)

Y_ij = matrix(0,N,N)
for(i in 1:N){
  for(j in 1:N){
    Y_ij[i,j] = rbinom(1,M, P_NbyN[i,j])
  }
}


z_k = which(z == 1)
P_NbyN[z_k,2]
mean(Y_ij[z_k,2]/M)


chains= list()
for(i in 1:4){
  alpha0=runif(1,0.1,3)
  trunc=improper_prior5(K,beta_max,alpha = alpha0)
  S0=runif(1,0.1,.9)
  P0_POMM= simulating_overlapping_POMM_powerlaw_norm(K,alpha0,S0,trunc,beta_max,diag0.5)
  init_pomm =list(z = rep_len(sample(1:K,K,F), N),alpha=alpha0,S=S,P=P0_POMM)
  names(init_pomm)
  estimation_control = list(z = 0,alpha=1,S=1,P=1)
  ground_truth= list(z = z,alpha=alpha,S=S,P=P)
  hyper_params = list(K = K,beta_max =beta_max,gamma_vec = gamma_vec,diag0.5=diag0.5)
  TEST = adaptive_MCMC_POMM(Yij_matrix = Y_ij,Nij_matrix = N_ij,init = init_pomm,estimation_control = estimation_control,ground_truth = ground_truth,N = N,N_iter = N_iter,targ_rate = .22,hyper_params =hyper_params ,seed = 123)
  chains[[i]]= TEST
}
#simple study

P0_Simple= matrix(.5,K,K)
P0_Simple[upper.tri(P0_Simple)]<- runif(K*(K-1)/2,0.5,beta_max)
P0_Simple[lower.tri(P0_Simple)]<- 1- P0_Simple[upper.tri(P0_Simple)]

init_Simple =list(z = rep_len(sample(1:K,K,F), N),P=P0_Simple)
names(init_Simple)
estimation_control_Simple = list(z = 1,P=1)
ground_truth_Simple= list(z = z,P=P)
hyper_params_Simple = list(K = K,beta_max =beta_max,gamma_vec = gamma_vec,diag0.5=diag0.5)
TEST = adaptive_MCMC_simple(Yij_matrix = Y_ij,Nij_matrix = N_ij,init = init_Simple,estimation_control = estimation_control_Simple,ground_truth = ground_truth_Simple,N = N,N_iter = N_iter,targ_rate = .22,hyper_params =hyper_params_Simple ,seed = 123)
chains[[1]]= TEST
}





mean(TEST$st.deviations$sd_p[1,2,])
mean(TEST$st.deviations$sd_p[1,3,])
mean(TEST$st.deviations$sd_p[2,3,])
for(i in 1:K){
  for(j in 1:K)
    print(mean(TEST$est_containers$P[i,j,]) - P[i,j])
}

mm<-mcmc.list(mcmc(TEST$est_containers$P[1,2,]),mcmc(TEST$est_containers$P[1,3,]),mcmc(TEST$est_containers$P[2,3,]))
acc <- list(TEST$acceptance_rates$acc.count_p,TEST$acceptance_rates$acc.count_p,TEST$acceptance_rates$acc.count_p)
gelman.diag(mm)

library(gt)
library(coda)

chains_list = list(TEST$est_containers$P[1,2,],TEST$est_containers$P[1,3,],TEST$est_containers$P[2,3,])
mcmc_fied<- (lapply(chains_list, mcmc))
mcmc_list_fied<- mcmc.list(mcmc_fied)

P_diagnostic_table(chains, T, diag0.5,P,K)

P_summary_table(test_output = test1 ,true_value =T,diag0.5 = T,K=K,P = P,burn_in = 3000)
z_summary_table(test1,model = 'POMM',burn_in = 2000)
z_

alpha_summary_table(test_output = test1,true_value = T,diag0.5 = T,alpha = alpha,K = K,burn_in = 3000)

P_diagnostic_table(chains = chains ,true_value = F,diag0.5 = T,K=K,P = P,burn_in = 3000)
z_diagnostic_table(chains = chains,true_value = F,diag0.5 = T,K=K,z = z,burn_in = 3000)
alpha_table<-alpha_diagnostic_table(chains = chains,true_value = F,diag0.5 = T,K=K,alpha=alpha,burn_in = 3000)

alpha_title<-paste0('alpha_',alpha,'K',K,'.tex')
alpha_table %>% gt()%>% as_latex()%>% as.character() %>% writeLines(con = alpha_title)
getwd()

test1$acceptance_rates$acc.count_alpha
# Define a function to save tables
save_table_to_file <- function(table_code, filename) {
  table_code %>%
    gt() %>%
    as_latex() %>%
    as.character() %>%
    writeLines(con = filename)
}

# Replace these values with your desired values for alpha and K
alpha <- 0.05
K <- 10

# Specify the directory where you want to save the files
# Replace "path/to/your/directory/" with the actual directory path
directory <- "path/to/your/directory/"

# Save each table with appropriate filenames and directory
alpha_title <- paste0(directory, 'alpha_', alpha, 'K', K, '.tex')

# Replace the table codes with the respective functions that generate the tables
P_summary_table(test_output = test1, true_value = FALSE, diag0.5 = TRUE, K = K, P = P, burn_in = 3000)
save_table_to_file(table_code, alpha_title)

z_summary_table(test1, model = 'POMM', burn_in = 2000)
save_table_to_file(table_code, paste0(directory, "z_summary_table.tex"))

alpha_summary_table(test_output = test4, true_value = TRUE, diag0.5 = TRUE, alpha = alpha, K = K, burn_in = 3000)
save_table_to_file(table_code, paste0(directory, "alpha_summary_table.tex"))

P_diagnostic_table(chains = chains, true_value = FALSE, diag0.5 = TRUE, K = K, P = P, burn_in = 3000)
save_table_to_file(table_code, paste0(directory, "P_diagnostic_table.tex"))

z_diagnostic_table(chains = chains, true_value = FALSE, diag0.5 = TRUE, K = K, z = z, burn_in = 3000)
save_table_to_file(table_code, paste0(directory, "z_diagnostic_table.tex"))

alpha_table <- alpha_diagnostic_table(chains = chains, true_value = FALSE, diag0.5 = TRUE, K = K, alpha = alpha, burn_in = 3000)
save_table_to_file(table_code, paste0(directory, "alpha_table.tex"))



View(S)

#chains[[i]]= TEST
#}

test1<-chains[[1]]
test2<-chains[[2]]
test3<-chains[[3]]
test4<-chains[[4]]


assembling_chains <- function(chains, burnin, parameter){
  test1<-chains[[1]]
  test2<-chains[[2]]
  test3<-chains[[3]]
  test4<-chains[[4]]
  if(parameter == 'S'){
    assembled<- c(test1$est_containers$S[-c(1:burnin)],test2$est_containers$S[-c(1:burnin)],test3$est_containers$S[-c(1:burnin)],test4$est_containers$S[-c(1:burnin)])
    return(assembled)
  }else if(parameter == 'alpha'){
    assembled<- c(test1$est_containers$alpha[-c(1:burnin)],test2$est_containers$alpha[-c(1:burnin)],test3$est_containers$alpha[-c(1:burnin)],test4$est_containers$alpha[-c(1:burnin)])
    return(assembled)
  }else if(parameter == 'P'){
    assembled<- abind::abind(test1$est_containers$P[,,-c(1:burnin)],test2$est_containers$P[,,-c(1:burnin)],test3$est_containers$P[,,-c(1:burnin)],test4$est_containers$P[,,-c(1:burnin)],along = 3)
    return(assembled)
  }else if(parameter == 'z'){
    assembled<- cbind(test1$est_containers$z[,-c(1:burnin)],test2$est_containers$z[,-c(1:burnin)],test3$est_containers$z[,-c(1:burnin)],test4$est_containers$z[,-c(1:burnin)])
    return(assembled)
  }else{
  print('Please provide a valide name for the parameters: alpha,S,P,Z')
  }
}


plot(TEST$st.deviations$sd_alpha)
plot(density(TEST$est_containers$alpha[-c(1:2000)]))

#extracting similarity matrix
similarity_matrixPOMMM = pr_cc(TEST$est_containers$z[,-c(1:1000)])
similarity_plot(similarity_matrixPOMMM, z, z) #checking mixing
library(coda)

#point est 1
point_est_POMM = minVI(similarity_matrixPOMMM,method = 'avg')$cl


S_container = assembling_chains(chains,burnin = 5000,parameter = 'S')
alpha_container = assembling_chains(chains,burnin = 5000,parameter = 'alpha')

data <- data.frame(alpha = alpha_container, S = S_container)
alpha_true=alpha
S_true<- S
# Create scatterplot with univariate distributions
ggplot(data, aes(x = alpha, y = S)) +
  geom_point(alpha=.8) +
  geom_point(x=alpha,y=S, col='red',size=2,shape=15)+
  geom_density_2d(size = .8, alpha=.7) +
  geom_xsidedensity(fill = 'red', alpha = 0.5) +
  geom_ysidedensity(fill = 'blue', alpha = 0.5) +
  labs(title = "Bivariate Scatterplot with Densities", x = "Alpha", y = "S") +
  theme_bw()





# Create scatterplot with univariate distributions
ggplot(data, aes(x = alpha, y = S)) +
  geom_point() +
  geom_xsidedensity(aes(fill = alpha), alpha = 0.5, position = "stack") +
  geom_ysidedensity(aes(fill = S), alpha = 0.5, position = "stack") +
  theme_minimal()



