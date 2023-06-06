


source("/Users/lapo_santi/Desktop/Nial/project/simplified model/Functions_priorSST.R")
improper_prior = function(K,alpha,beta_max){
  x = c(1:K)
  r = x^alpha
  r_scaled=(r-min(r))/(max(r)-min(r))*(beta_max-0.5)+0.5
  return(r_scaled)
}

improper_prior2 = function(K,alpha,beta_max){
  x = c(1:K)
  r = x^alpha
  r_scaled=(r-min(x))/(max(x)-min(x))*(beta_max-0.5)+0.5
  r = r_scaled^alpha
  return(r_scaled)
}

improper_prior3 = function(K,beta_max,alpha){
  x_min = 0.5^(1/alpha)
  x_max = beta_max^(1/alpha)
  delta = (x_max - x_min)/K
  delta_sum = rep(delta,K-1)
  points = cumsum(c(x_min,delta_sum))
  boundaries = points^alpha
  return(boundaries)}

improper_prior4 = function(K,beta_max,alpha){
  x_min =  exp(.5/alpha) - 1
  x_max = exp(beta_max/alpha) - 1
  delta = (x_max - x_min)/K
  delta_sum = rep(delta,K-1)
  points = cumsum(c(x_min,delta_sum))
  boundaries = alpha*(log(points+1))
  return(boundaries)}

beta_max = sample_beta_trunc(1,1,1,0.5,1,6)
alpha = sample_norm_trunc(1,1,1,0,Inf)







################################################
#Generating a matrix according to this new model
###############################################


###############################################
#Plotting the heatmap of a single configuration
###############################################





beta_max = sample_beta_trunc(1,1,1,0.5,1,6)
alpha = sample_norm_trunc(1,1,1,0,Inf)


K=10
alpha=3
beta_max=.8

###############
# plot of the boundary process
###############
alpha_simulation = c(1/10, 1/5,1/3,1/2,1,1.5,2,3)

bounds_container = matrix(0,nrow = length(alpha_simulation), ncol=K+1)
for( i in 1:length(alpha_simulation)){
  sim_i = improper_prior5(K,beta_max = .9,alpha_simulation[i])
  bounds_container[i,] = sim_i
}


# Plot the lines
ts.plot(bounds_container[1,], col = "red", xlab = "Level Sets", ylab = "Boundaries", main = "Boundaries evolution process")
for(i in 2:nrow(bounds_container)) {
  lines(bounds_container[i,], col = rainbow(length(alpha_simulation))[i], type = "l")
}

# Add a legend
legend("bottomright", legend = round(alpha_simulation,2), col = rainbow(length(alpha_simulation)), lty = 1, cex=.45)

###############
# plot of the P matrix 
###############


K=10
alpha=1.5
beta_max=.8
test = simulating_POMM_powerlaw(K,alpha,beta_max)

heat_map_blue = function(matrix, title){
# convert matrix to data frame
df <- as.data.frame(t(matrix))
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
  labs(title=title, x="Column", y="Row", fill="Prob")+
  geom_text(aes(label = round(Z, 2)), size =2)
}


# Load ggplot2 library
library(ggplot2)

# Define alpha_simulation and bounds_container
alpha_simulation = c(1/10, 1/5, 1/3, 1/2, 1, 1.5, 2, 3)
bounds_container = matrix(0, nrow = length(alpha_simulation), ncol = K+1)

# Calculate bounds_container for each alpha_simulation
for(i in 1:length(alpha_simulation)) {
  sim_i = improper_prior5(K, beta_max = 0.9, alpha_simulation[i])
  bounds_container[i, ] = sim_i
}

# Convert bounds_container to a data frame for ggplot2
df <- data.frame(Time = c(1:11), t(bounds_container))
colnames(df) <-  c("Time","1/10" ,"1/5", "1/3","1/2",'1','2',"3/2","3")

# Create the plot using ggplot2
ggplot(df, aes(x = Time)) +
  geom_line(aes(y = `1/10`, color = "1/10")) +
  geom_line(aes(y = `1/5`, color = "1/5")) +
  geom_line(aes(y = `1/3`, color = "1/3")) +
  geom_line(aes(y = `1/2`, color = "1/2")) +
  geom_line(aes(y = `1`, color = "1")) +
  geom_line(aes(y = `2`, color = "3/2")) +
  geom_line(aes(y = `3/2`, color = "2")) +
  geom_line(aes(y = `2`, color = "3")) +
  xlab("Level Sets") +
  ylab("Boundaries") +
  ggtitle("Boundaries evolution") +
  scale_color_manual(values = rainbow(length(alpha_simulation))) +
  theme(legend.position = "bottomright", legend.title = element_blank(), legend.text = element_text(size = 8))



pdf_power_law <- function(y, k, K, beta_max) {
  abs(2/(log((k-1)/K)*(2*y-1))) *
    (1/sqrt(2*pi)) *
    exp(-0.5*( (log(y-0.5) - log(beta_max-0.5))/(log((k-1)/K)) - 1.5)^2)
}
try = improper_prior5(10,.9,1.5)









k = c(2:10)
for(i in k ){
x= seq(.51,1,.01)
y=pdf_power_law(x,k,10,.9)   
plot(x,y,"l")
abline(v= try[k])
}
# Define the range of x values
x <- seq(.5, 1, .01)

# Set up the plotting window
par(mfrow = c(3, 3), mar = c(3, 3, 1, 1), oma = c(0, 0, 2, 0))

# Loop over different k values
for (i in 2:K) {
  # Calculate the PDF for the current k value
  y <- pdf_power_law(x, i, 10, .9)
  
  # Plot the PDF
  plot(x, y, type = "l", xlab = "x", ylab = "PDF", main = paste0("i = ", i))
  
  # Add a vertical line at the current k value
  abline(v = try[i], col = "red")
}

# Add a main title to the entire page
mtext("PDFs for different each level set i, alpha=1.5", outer = TRUE, line = 1)


integrate(pdf_power_law,.5,10, K=10, beta_max=.9, k=4)

