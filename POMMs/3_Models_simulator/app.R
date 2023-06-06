#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(plotly)

source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")

semi_symmetric = function(Y){
  Y = (1 -t(Y*upper.tri(Y))) * (lower.tri(Y, diag = F)*1) + upper.tri(Y, diag = T) *Y
  return(Y)}



plot_matrix <- function(my_matrix, cilower, ciupper, model_name, colorscale, width = 400, height = 300) {
  n <- nrow(my_matrix)
  plot_ly(x = 1:n, y = 1:n, z = my_matrix, type = "heatmap", colorscale = colorscale, width = width, height = height) %>%
    #add_surface(x = 1:n, y = 1:n, z = cilower, opacity = 0.4, showscale = FALSE) %>%
    #add_surface(x = 1:n, y = 1:n, z = ciupper, opacity = 0.4, showscale = FALSE) %>%
    colorbar(title = "Mean") %>%
    layout(title = model_name,
           xaxis = list(title = "Column index"),
           yaxis = list(title = "Row index"))
}

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



# Define UI for app
ui <- shinyUI(fluidPage(
  titlePanel("Simulation of three different POMMs models "),
  sidebarLayout(
    sidebarPanel(
      sliderInput("K", "Number of teams/players:", min = 0, max = 100, value = 10, step = 1),
      sliderInput("lambda_0", "Lambda:", min = 0.001, max = 0.05, value = 0.01, step = 0.001),
      sliderInput("N_iter", "Number of Iterations:", min = 0, max = 10000, value = 500, step = 1),
      actionButton("runButton", "Run")
    ),
    mainPanel(
      "Main Panel",
      fluidRow(
        splitLayout(cellWidths = c("33%", "33%","33%"), 
                    plotlyOutput("plotgraph1"), 
                    plotlyOutput("plotgraph2"), 
                    plotlyOutput("plotgraph3"))
      )
    )
  )
))

# Define server logic required to draw plots
server <- function(input, output, session) {
  source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R", local = T)
  
  # Render plot 1
  output$plotgraph1 <- renderPlotly({
    # Set up the parameters
    n <- input$K
    lambda_0 <- input$lambda_0
    N_iter <- input$N_iter
  
    # Initialize the result matrices
    means1 <- matrix(NA, n, n)
    cilower1 <- matrix(NA, n, n)
    ciupper1 <- matrix(NA, n, n)
    
    
    
    values1 <- replicate(N_iter,simulating_POMM_model_1(n, lambda_0 = lambda_0))
    for (i in 1:n) {
      for (j in 1:i) {
        means1[i, j] <- mean(values1[i,j,],na.rm = T)
        means1[j, i] <- mean(values1[j,i,],na.rm = T)
        cilower1[i, j] <- quantile(values1[i,j,], 0.025,na.rm = T)
        cilower1[j, i] <- quantile(values1[j,i,], 0.025,na.rm = T)
        ciupper1[i, j] <- quantile(values1[i,j,], 0.975,na.rm = T)
        ciupper1[j, i] <- quantile(values1[j,i,], 0.975,na.rm = T)}}
    
    
    plot_matrix(means1, cilower1, ciupper1, "Model 1", "Viridis")
    })
  
  # Render plot 2
  output$plotgraph2 <- renderPlotly({
    # Set up the parameters
    n <- input$K
    lambda_0 <- input$lambda_0
    N_iter <- input$N_iter
    
    
    # Initialize the result matrices
    
    means2 <- matrix(NA, n, n)
    cilower2 <- matrix(NA, n, n)
    ciupper2 <- matrix(NA, n, n)
    
    values2 <- replicate(N_iter,simulating_POMM_model_2(n, lambda_0 = lambda_0))
    for (i in 1:n) {
      for (j in 1:i) {
        means2[i, j] <- mean(values2[i,j,],na.rm = T)
        means2[j, i] <- mean(values2[j,i,],na.rm = T)
        cilower2[i, j] <- quantile(values2[i,j,], 0.025,na.rm = T)
        cilower2[j, i] <- quantile(values2[j,i,], 0.025,na.rm = T)
        ciupper2[i, j] <- quantile(values2[i,j,], 0.975,na.rm = T)
        ciupper2[j, i] <- quantile(values2[j,i,], 0.975,na.rm = T)}}
    
    
       plot_matrix(means2, cilower2, ciupper2, "Model 2", "Viridis") })
  # Render plot 3
  output$plotgraph3 <- renderPlotly({
    # Set up the parameters
    n <- input$K
    lambda_0 <- input$lambda_0
    N_iter <- input$N_iter
    
    
    # Initialize the result matrices
    
    means3 <- matrix(NA, n, n)
    cilower3<- matrix(NA, n, n)
    ciupper3<- matrix(NA, n, n)
    
    
    values3 <- replicate(N_iter,simulating_POMM_model_3(n, lambda_0 = lambda_0))
    for (i in 1:n) {
      for (j in 1:i) {
        means3[i, j] <- mean(values3[i,j,],na.rm = T)
        means3[j, i] <- mean(values3[j,i,],na.rm = T)
        cilower3[i, j] <- quantile(values3[i,j,], 0.025,na.rm = T)
        cilower3[j, i] <- quantile(values3[j,i,], 0.025,na.rm = T)
        ciupper3[i, j] <- quantile(values3[i,j,], 0.975,na.rm = T)
        ciupper3[j, i] <- quantile(values3[j,i,], 0.975,na.rm = T)}}
    
    
    plot_matrix(means3, cilower3, ciupper3, "Model 3", "Viridis") })
  
}

shinyApp(ui = ui, server = server)
  
  
  