library(phyloseq)
library(vegan)
library(ggplot2)
library(reshape2)

"""
This script processes a single sample within a phyloseq analysis workflow. 
It performs the following tasks:

1. Loads species data from a provided CSV file.
2. Normalizes the species data to relative abundances.
3. Generates a CSV file with the normalized species counts.
4. Creates and saves a richness plot, which visualizes the number of species (richness) in the sample.

The script takes three command-line arguments:
    1. data_file: The path to the input CSV file containing species data.
    2. richness_plot: The path where the richness plot will be saved.
    3. species_count_file: The path where the normalized species count file will be saved.
"""

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
diversity_plot <- args[2]
species_count_file <- args[3]

otu_data <- read.csv(data_file, row.names = 1)

if (nrow(otu_data) == 0) {
    quit(status = 0)
}

otu_data_normalized <- otu_data / rowSums(otu_data)
write.csv(otu_data_normalized, file = species_count_file)

# Calculate diversity metrics
observed_species <- sum(rowSums(otu_data) > 0)  # OTUs count
chao1_index <- estimateR(otu_data)["S.chao1", 1]  # Chao1 index
fisher_alpha <- fisher.alpha(rowSums(otu_data))  # Fisher-Alpha
shannon_index <- diversity(otu_data, index = "shannon")  # Shannon diversity
simpson_index <- 1 - diversity(otu_data, index = "simpson")  # Simpson diversity, converted to 1 - Simpson

# Combine metrics into a data frame
diversity_metrics <- data.frame(
    Metric = c("OTUs", "Chao1 Index", "Fisher-Alpha", "Shannon Index", "Simpson Index"),
    Value = c(observed_species, chao1_index, fisher_alpha, shannon_index, simpson_index)
)

# Melt the data for ggplot2
diversity_metrics_melt <- melt(diversity_metrics, id.vars = "Metric")

# Plot using ggplot2
png(diversity_plot, width = 1200, height = 400)  # Adjust size to accommodate horizontal layout
ggplot(diversity_metrics_melt, aes(x = Metric, y = value, fill = Metric)) +
    geom_bar(stat = "identity", color = "black") +
    geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.5, size = 3) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    ylab("Index Value") +
    ggtitle("Diversity Indices") +
    facet_wrap(. ~ Metric, scales = "free_y", nrow = 1) +  # Arrange all plots in a single row
    theme(axis.text.x = element_blank(),  # Remove x-axis text
          axis.ticks.x = element_blank())  # Remove x-axis ticks
dev.off()
