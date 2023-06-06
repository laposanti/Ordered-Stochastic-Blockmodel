library(MASS)

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu^2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

diag_matrix = function(K){
  my_diagonal_matrix  = matrix(NA,K,K)
  my_diagonal_matrix = row(my_diagonal_matrix) - col(my_diagonal_matrix)
  return(my_diagonal_matrix)
}

diag_split_NA = function(K){
  aux_matrix= matrix(NA,K,K)
  diag_indicator <- (row(aux_matrix) - col(aux_matrix))
  diag_indicator<- upper.tri(diag_indicator,diag = T)*diag_indicator + matrix(-1,nrow(diag_indicator),ncol(diag_indicator))*upper.tri(diag_indicator,diag = T)
  L_k = split(aux_matrix, diag_indicator)
  return(L_k)
}
diag_split_matrix = function(matrix_0){
  K = nrow(matrix_0)
  diag_indicator <- (row(matrix_0) - col(matrix_0))
  diag_indicator<- upper.tri(diag_indicator,diag = T)*diag_indicator + matrix(-1,nrow(diag_indicator),ncol(diag_indicator))*upper.tri(diag_indicator,diag = T)
  L_k = split(matrix_0, diag_indicator)
  return(L_k)
}

sample_beta_trunc = function(N,alpha,beta,a,b){
  u = runif(N)
  return( qbeta(((pbeta(b, alpha, beta) -pbeta(a, alpha, beta))*u+pbeta(a, alpha, beta)),alpha,beta))}

K=4
simulating_POMM_model_1 = function(K, lambda_0,  alpha_0 = .5){
  #simulating from the prior
  #lambda_0 = rgamma(1,shape = alpha_0, scale = lambda_0)

  #L_k is a list containing each diagonal of the upper triangular matrix
  L_k = diag_split_NA(K)
  L_k[[K]] = rep(0.5,K)
 
  # Compute beta values
  for(i in (K-1):1){
    n_items = length(L_k[[i]])
    max_adj = matrix(NA, nrow=n_items,ncol= 1)
    for(j in 1:n_items){
      #obtaining the maximum of the 2 adjl neighborhoods
      max_adj[j] = max(L_k[[i+1]][j], L_k[[i+1]][j+1])
    }
    #max of the adjil + the gamma advantage factor
    mu = max_adj 
    #the variance is constant
    
    var= min(lambda_0, (max_adj*(1-max_adj)-0.0001))
    #setting alpha and beta s.t. E(X) = M_adjl + gamma_ij
    alpha= estBetaParams(mu, var)$alpha
    beta = estBetaParams(mu, var)$beta
    #generating the beta values
    L_k[[i]] = sample_beta_trunc(n_items,alpha,beta,max_adj, .99)
  }
  #putting back everything together in the matrix
  result = diag_matrix(K)
  for(i in 0:(K-1)){
    result[which(result==-i)]<-L_k[[K-i]]
  }
  result = semi_symmetric(result)
  return(result)}

alpha_0 = .5
lambda_0=.2
rgamma(1,shape = alpha_0, scale = lambda_0)
simulating_POMM_model_3(4, lambda_0=.04)

simulating_POMM_model_2 = function(K, lambda_0,  alpha_0 = .5){
  #simulating from the prior
  #lambda_0 = rgamma(1,shape = alpha_0, scale = lambda_0)
  
  #L_k is a list containing each diagonal of the upper triangular matrix
  L_k = diag_split_NA(K)
  L_k[[K]] = rep(0.5,K)
  
  
  #diagonal matrix for proportionality
  
  prop_matrix <- diag_split_matrix(diag_matrix(K))
  # Compute beta values
  for(i in (K-1):1){
    n_items = length(L_k[[i]])
    max_adj = matrix(NA, nrow=n_items,ncol= 1)
    for(j in 1:n_items){
      #obtaining the maximum of the 2 adjl neighborhoods
      max_adj[j] = max(L_k[[i+1]][j], L_k[[i+1]][j+1])
    }
    #max of the adjil + the gamma advantage factor
    mu = max_adj 
    #the variance is constant

    var= min(lambda_0*(-prop_matrix[[i]]), (max_adj*(1-max_adj)-0.0001))
    #setting alpha and beta s.t. E(X) = M_adjl + gamma_ij
    alpha= estBetaParams(mu, var)$alpha
    beta = estBetaParams(mu, var)$beta
    #generating the beta values
    L_k[[i]] = sample_beta_trunc(n_items,alpha,beta,max_adj, .99)
  }
  #putting back everything together in the matrix
  result = diag_matrix(K)
  for(i in 0:(K-1)){
    result[which(result==-i)]<-L_k[[K-i]]
  }
  result = semi_symmetric(result)
  return(result)}

simulating_POMM_model_3 = function(K, lambda_0,  alpha_0 = .5){
  #simulating from the prior
  gamma_0 = rgamma(K-1, shape = alpha_0, scale = lambda_0)
  #lambda_0 = rgamma(1,shape = alpha_0, scale = lambda_0)
  
  #L_k is a list containing each diagonal of the upper triangular matrix
  L_k = diag_split_NA(K)
  L_k[[K]] = rep(0.5,K)
  
  
  #diagonal matrix for proportionality
  prop_matrix <- diag_split_matrix(diag_matrix(K))
  # Compute beta values
  for(i in (K-1):1){
    n_items = length(L_k[[i]])
    max_adj = matrix(NA, nrow=n_items,ncol= 1)
    for(j in 1:n_items){
      #obtaining the maximum of the 2 adjl neighborhoods
      max_adj[j] = max(L_k[[i+1]][j], L_k[[i+1]][j+1])
    }
    #max of the adjil + the gamma advantage factor
    mu = max_adj 
    #the variance is constant
    
    var= min(gamma_0[i]*(-prop_matrix[[i]]), (max_adj*(1-max_adj)-0.0001))
    #setting alpha and beta s.t. E(X) = M_adjl + gamma_ij
    alpha= estBetaParams(mu, var)$alpha
    beta = estBetaParams(mu, var)$beta
    #generating the beta values
    L_k[[i]] = sample_beta_trunc(n_items,alpha,beta,max_adj, .99)
  }
  #putting back everything together in the matrix
  result = diag_matrix(K)
  for(i in 0:(K-1)){
    result[which(result==-i)]<-L_k[[K-i]]
  }
  result = semi_symmetric(result)
  return(result)}




# Generate matrices for each model
set.seed(123)
n=4
# Calculate means and confidence intervals for each element in the matrices
n <- nrow(model1)
means1 <- matrix(NA, n, n)
means2 <- matrix(NA, n, n)
means3 <- matrix(NA, n, n)
cilower1 <- matrix(NA, n, n)
ciupper1 <- matrix(NA, n, n)
cilower2 <- matrix(NA, n, n)
ciupper2 <- matrix(NA, n, n)
cilower3 <- matrix(NA, n, n)
ciupper3 <- matrix(NA, n, n)

N_iter <- 10000
values1 <- replicate(N_iter,simulating_POMM_model_1(n, lambda_0 = 0.02))
values2 <- replicate(N_iter,simulating_POMM_model_2(n, lambda_0 = 0.02))
values3 <- replicate(N_iter,simulating_POMM_model_3(n, lambda_0 = 0.02))
for (i in 1:n) {
  for (j in 1:i) {
    means1[i, j] <- mean(values1[i,j,],na.rm = T)
    means1[j, i] <- mean(values1[j,i,],na.rm = T)
    cilower1[i, j] <- quantile(values1[i,j,], 0.025,na.rm = T)
    cilower1[j, i] <- quantile(values1[j,i,], 0.025,na.rm = T)
    ciupper1[i, j] <- quantile(values1[i,j,], 0.975,na.rm = T)
    ciupper1[j, i] <- quantile(values1[j,i,], 0.975,na.rm = T)
    
    means2[i, j] <- mean(values2[i,j,],na.rm = T)
    means2[j, i] <- mean(values2[j,i,],na.rm = T)
    cilower2[i, j] <- quantile(values2[i,j,], 0.025,na.rm = T)
    cilower2[j, i] <- quantile(values2[j,i,], 0.025,na.rm = T)
    ciupper2[i, j] <- quantile(values2[i,j,], 0.975,na.rm = T)
    ciupper2[j, i] <- quantile(values2[j,i,], 0.975,na.rm = T)
    
    means3[i, j] <- mean(values3[i,j,],na.rm = T)
    means3[j, i] <- mean(values3[j,i,],na.rm = T)
    cilower3[i, j] <- quantile(values3[i,j,], 0.025,na.rm = T)
    cilower3[j, i] <- quantile(values3[j,i,], 0.025,na.rm = T)
    ciupper3[i, j] <- quantile(values3[i,j,], 0.975,na.rm = T)
    ciupper3[j, i] <- quantile(values3[j,i,], 0.975,na.rm = T)
  }
}

from_matrix_to_df = function(my_matrix){
  K= nrow(my_matrix)
  means1_df = vector()
  for( i in 1:(K-1)){
    means1_df = append(means1_df,split(my_matrix,diag_matrix(K))[[i]])
  }
  return(means1_df)
}

df_model1 = data.frame(mean =from_matrix_to_df(means1), ci_lower = from_matrix_to_df(cilower1),ci_upper = from_matrix_to_df(ciupper1))
df_model2 = data.frame(mean =from_matrix_to_df(means2), ci_lower = from_matrix_to_df(cilower2),ci_upper = from_matrix_to_df(ciupper2))
df_model3 = data.frame(mean =from_matrix_to_df(means3), ci_lower = from_matrix_to_df(cilower3),ci_upper = from_matrix_to_df(ciupper3))





plot_matrix <- function(matrix, cilower, ciupper, model_name, colorscale) {
  n <- nrow(matrix)
  plot_ly(x = 1:n, y = 1:n, z = matrix, type = "surface", colorscale = colorscale) %>%
    add_surface(x = 1:n, y = 1:n, z = cilower, opacity = 0.4, showscale = FALSE) %>%
    add_surface(x = 1:n, y = 1:n, z = ciupper, opacity = 0.4, showscale = FALSE) %>%
    colorbar(title = "Mean") %>%
    layout(title = model_name,
           xaxis = list(title = "Column index"),
           yaxis = list(title = "Row index"))
}

plot1 <- plot_matrix(means1, cilower1, ciupper1, "1", "Viridis") %>% layout(layout)
plot2 <- plot_matrix(means2, cilower2, ciupper2, "2", "Brewer Blues")
plot3 <- plot_matrix(means3, cilower3, ciupper3, "3", "Viridis")
layout <- list(
  title = 'Side By Side Subplots',
  showlegend = TRUE,
  scene = list(
    xaxis = list(title = 'Column index'),
    yaxis = list(title = 'Row index'),
    zaxis = list(title = 'Value')
  )
)
fig <- subplot(plot1, plot2, plot3)  %>% layout(layout)

fig

plot1
plot2
plot3
#######
#Not working code
#######
# library(plotly)
# 
# # Define the three plots
# plot1 <- plot_matrix(means1, cilower1, ciupper1, "1")
# plot2 <- plot_matrix(means2, cilower2, ciupper2, "2")
# plot3 <- plot_matrix(means3, cilower3, ciupper3, "3")
# 
# # Define the layout of the subplots
# layout <- list(
#   title = 'Side By Side Subplots',
#   showlegend = FALSE,
#   annotations = list(
#     list(x = 0.16, y = 1, xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom", text = "1"),
#     list(x = 0.5, y = 1, xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom", text = "2"),
#     list(x = 0.84, y = 1, xref = "paper", yref = "paper", xanchor = "center", yanchor = "bottom", text = "3")
#   ),
#   xaxis = list(domain = c(0, 0.33), showgrid = FALSE),
#   yaxis = list(domain = c(0.5, 1), showgrid = FALSE),
#   xaxis2 = list(domain = c(0.34, 0.66), showgrid = FALSE),
#   yaxis2 = list(domain = c(0.5, 1), showgrid = FALSE),
#   xaxis3 = list(domain = c(0.67, 1), showgrid = FALSE),
#   yaxis3 = list(domain = c(0.5, 1), showgrid = FALSE)
# )



library(ggplot2)
library(cowplot)
plot_matrix <- function(matrix, cilower, ciupper, model_name) {
  n <- nrow(matrix)
  data <- expand.grid(x = 1:n, y = 1:n)
  data$z <- matrix[cbind(data$x, data$y)]
  data$cilower <- cilower[cbind(data$x, data$y)]
  data$ciupper <- ciupper[cbind(data$x, data$y)]
  
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_tile(aes(fill = z)) +
    scale_fill_gradient(name = "Mean") +
    labs(title = model_name, x = "Column index", y = "Row index") +
    theme_bw()
  
  p <- p + geom_tile(aes(fill = mean), alpha = 0.5) + 
    guides(fill = guide_legend(title = "CI", override.aes = list(alpha = 1))) +
    geom_tile(aes(fill = ciupper), alpha = 0.5)
  
  p
}

plot1 <- plot_matrix(means1, cilower1, ciupper1, "1")
plot2 <- plot_matrix(means2, cilower2, ciupper2, "2")
plot3 <- plot_matrix(means3, cilower3, ciupper3, "3")

fig <- plot_grid(plot1, plot2, plot3, ncol = 3)
fig


