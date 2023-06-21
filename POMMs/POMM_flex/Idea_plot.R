# Install required packages
install.packages("scatterplot3d")


# Load required packages
library(scatterplot3d)
library(rgl)
# Create three normal distributions
set.seed(123)  # Set seed for reproducibility
x1 <- rnorm(1000, mean = 0.5, sd = .1)
x2 <- rnorm(1000, mean = 0.6, sd = .1)
x3 <- rnorm(1000, mean = 0.7, sd = .1)

# Compute the density values for each x value
d1 <- density(x1)
d2 <- density(x2)
d3 <- density(x3)

# Create a data frame with x, y, and z values
df <- data.frame(
  x = c(d1$x, d2$x, d3$x),
  y = rep(1:3, each = length(d1$x)),
  z = c(d1$y, d2$y, d3$y)
)

# Create a blue color scale
color_scale <- colorRampPalette(c("lightblue", "darkblue"))(3)
color_scale<- c(rep(color_scale[1], length(d1$x)), rep(color_scale[2], length(d2$x)),rep(color_scale[3], length(d3$x)))
# Create the 3D plot with rotation and color scale
scatterplot3d(
  df$x, df$y, df$z,
  xlab = "X Values",
  ylab = "Discrete Scale",
  xlim = c(0,1),
  zlab = "Density",
  pch = 16,
  type = "h",
  box=F,
  color = color_scale, 
  angle=210,
  lwd = 0.4,      # Adjust line width
)

