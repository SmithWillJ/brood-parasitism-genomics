##################################
### Plotting admixture results ###
##################################

# Set working directory

setwd("/filepath/ReeWa_Admix")

# Loading required libraries

library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

# Read info files

coords <- read.table("../Individuals.coords", header = FALSE, col.names = c("lon", "lat"))
popmap <- read.table("../PopMap", header = FALSE, col.names = c("IndName", "Location"))

# Add individual index to coords (row order matches Q file order)

coords$IndName <- popmap$IndName
coords$Location <- popmap$Location

# Read and process Q files

read_q_file <- function(k_value) {
  filepath <- paste0("admixture_K", k_value, "/GF_auto_0.5.admix.", k_value, ".Q")
  data <- read.table(filepath, header = FALSE)
  
  # Add sample names from popmap
  data$IndName <- popmap$IndName
  data$K <- k_value
  
  return(data)
}

# Read Q files for K2, K3, and K4

k2_data <- read_q_file(2)
k3_data <- read_q_file(3)
k4_data <- read_q_file(4)

# Prepare data for plotting

prepare_plot_data <- function(data, k_value, coords_data) {
 
   # Reshape data to long format
  
  data_long <- data %>%
    select(-K) %>%
    pivot_longer(cols = -IndName, 
                 names_to = "Ancestry", 
                 values_to = "Proportion") %>%
    mutate(K = k_value,
           Ancestry = paste0("Ancestry_", as.numeric(gsub("V", "", Ancestry)))) %>%
   
     # Join with coordinates and location info
    
    left_join(coords_data %>% select(IndName, Location, lat), by = "IndName")
  
  return(data_long)
}

# Prepare data

k2_plot_data <- prepare_plot_data(k2_data, 2, coords)
k3_plot_data <- prepare_plot_data(k3_data, 3, coords)
k4_plot_data <- prepare_plot_data(k4_data, 4, coords)

# Create ordered individual factor (by latitude, from north to south)

individual_order_df <- coords %>%
  arrange(desc(lat)) %>%
  select(IndName, Location, lat)

individual_order <- individual_order_df$IndName

# Print the ordered list to console

cat("Individual Order (North to South by Latitude):\n")
cat("==============================================\n")
print(individual_order_df)
cat("\n")

# Save the ordered list to a file

write.csv(individual_order_df, "individual_order_by_latitude.csv", row.names = FALSE)
cat("Order saved to: individual_order_by_latitude.csv\n\n")

# Convert individual to factor with proper ordering

k2_plot_data$IndName <- factor(k2_plot_data$IndName, levels = individual_order)
k3_plot_data$IndName <- factor(k3_plot_data$IndName, levels = individual_order)
k4_plot_data$IndName <- factor(k4_plot_data$IndName, levels = individual_order)

# Define colours for ancestry components

ancestry_colours <- c(
  "Ancestry_1" = "#E41A1C",
  "Ancestry_2" = "#377EB8",
  "Ancestry_3" = "#4DAF4A",
  "Ancestry_4" = "#984EA3"
)

# Create individual plots with small gaps between bars

plot_k2 <- ggplot(k2_plot_data, aes(x = IndName, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity", width = 0.95) +
  scale_fill_manual(values = ancestry_colours) +
  labs(title = "K = 2", x = "", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right")

plot_k3 <- ggplot(k3_plot_data, aes(x = IndName, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity", width = 0.95) +
  scale_fill_manual(values = ancestry_colours) +
  labs(title = "K = 3", x = "", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right")

plot_k4 <- ggplot(k4_plot_data, aes(x = IndName, y = Proportion, fill = Ancestry)) +
  geom_bar(stat = "identity", width = 0.95) +
  scale_fill_manual(values = ancestry_colours) +
  labs(title = "K = 4", x = "Sample", y = "Ancestry Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right")

# Combine plots vertically

combined_plot <- gridExtra::grid.arrange(plot_k2, plot_k3, plot_k4, ncol = 1)

# Save the plot

png("admixture_K2_K3_K4.png", width = 1200, height = 900, res = 100)
print(combined_plot)
dev.off()
