# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
processed_wd<-"/Users/lapo_santi/Desktop/Nial/POMM_pairwise/POMMs/results/MCMC_output/Fixed_K/Application/Tennis_processed/"
# Sample data: Replace this with your actual data
set.seed(123)
Y_ij <- readRDS("./Data/Tennis application/Y_new.rds")
entities <- rownames(Y_ij)
z_SST <- read.csv(paste0(processed_wd,"z_est_K_6modelSST.csv"))[,2]
z_Simple <-  read.csv(paste0(processed_wd,"z_est_K_6modelSimple.csv"))[,2]

# Create a dataframe
df <- data.frame(entity = entities, z_SST = z_SST, z_Simple = z_Simple)

# Create y-axis positions based on z_est1 and z_est2
spacing_factor <- 1.5
df <- df %>%
  mutate(marginal_wins = rowSums(Y_ij)/(colSums(Y_ij)+rowSums(Y_ij)))%>%
  mutate(reorder_zSST = z_SST-marginal_wins)%>%
  mutate(reorder_zSimple = z_Simple-marginal_wins)%>%
  arrange(reorder_zSST, entity) %>%
  mutate(ySST = row_number()*spacing_factor) %>% # y-axis for z_est1
  arrange(reorder_zSimple, entity) %>%
  mutate(ySimple = row_number()*spacing_factor)     # y-axis for z_est2

# Transform the dataframe into a long format
df_long <- df %>%
  pivot_longer(cols = starts_with("z_"),
               names_to = "Partition",
               values_to = "Cluster") %>%
  mutate(Partition = factor(Partition, levels = c("z_SST", "z_Simple")))

# Plot using ggplot
ggplot(df_long, aes(x = Partition, y = reorder(entity,Cluster), group = entity)) + # Line connecting the same entity across partitions
  geom_point(aes(color = factor(Cluster)), size = 4) +   # Points colored by cluster            # Ensure entities are labeled correctly on y-axis
  labs(x = "Partition", y = "Entity", color = "Cluster") +
  theme_minimal() +
  theme(legend.position = "bottom") 




# Plot using ggplot
ggplot(df) +
  geom_segment(aes(x = 1, xend = 2, y = ySST, yend = ySimple, color = factor(z_SST)), size = 1, alpha=0.6,) + # Connecting lines
  geom_point(aes(x = 1, y = ySST, color = factor(z_SST)), size = 3) +                       # Points for z_est1
  geom_point(aes(x = 2, y = ySimple, color = factor(z_Simple)), size = 3) +                       # Points for z_est2
  scale_y_continuous(breaks = NULL) +                                                     # Remove default y-axis labels
  scale_x_continuous(breaks = c(1, 2), labels = c("Partition SST", "Partition Simple"),limits =c(-1,5) ) +        # Label the partitions
  labs(x = "", y = "", color = "Cluster") +
  theme_minimal() +
  geom_text(aes(x = 0.95, y = ySST, label = entity), hjust = 1) +                            # Labels for left y-axis
  geom_text(aes(x = 2.05, y = ySimple, label = entity), hjust = 0)+
  theme(legend.position.inside= c(1.5,1))

# Labels for right y-axis















