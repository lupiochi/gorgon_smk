library(vegan)
library(ggplot2)

"""
This script processes a single sample within a diversity analysis workflow. 
It performs the following tasks:

1. Loads species data from a provided CSV file.
2. Normalizes the species data to relative abundances.
3. Generates a CSV file with the normalized species counts.
4. Calculates diversity indices (Shannon and Simpson).
5. Creates and saves a diversity plot, which visualizes the Shannon and Simpson indices.

The script takes three command-line arguments:
    1. data_file: The path to the input CSV file containing species data.
    2. diversity_plot: The path where the diversity plot will be saved.
    3. species_count_file: The path where the normalized species count file will be saved.
"""

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
diversity_plot <- args[2]
species_count_file <- args[3]

# Read and process data
otu_data <- read.csv(data_file, row.names = 1)

# Check if the data is empty and exit if so
if (nrow(otu_data) == 0) {
    quit(status = 0)
}

# Normalize the data
otu_data_normalized <- otu_data / rowSums(otu_data)
write.csv(otu_data_normalized, file = species_count_file)

# Calculate diversity metrics
shannon <- mean(diversity(otu_data, index = "shannon"))
simpson <- mean(1 - diversity(otu_data, index = "simpson"))

# Create data frame for plotting
plot_data <- data.frame(
  Metric = c("Shannon", "Simpson"),
  Value = c(shannon, simpson)
)

# Create the plot
p <- ggplot(plot_data, aes(x = Metric, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", Value)), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Shannon" = "#ff7f0e", "Simpson" = "#2ca02c")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_blank()
  ) +
  labs(title = "Diversity Indices",
       x = NULL,
       y = "Index Value") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 1))

# Save the plot
ggsave(diversity_plot, plot = p, width = 6, height = 5, dpi = 300, bg = "white")

# Print the values
print(plot_data)