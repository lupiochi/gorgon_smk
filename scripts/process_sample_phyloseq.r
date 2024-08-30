library(phyloseq)

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
if (length(args) < 3) {
    stop("Insufficient arguments provided. Usage: Rscript process_sample_phyloseq.r <data_file> <richness_plot> <species_count_file>")
}
data_file <- args[1]
richness_plot <- args[2]
species_count_file <- args[3]

# Load the single sample data
otu_data <- read.csv(data_file, row.names = 1)

# Ensure the data has at least one row (species)
if (nrow(otu_data) == 0) {
    cat("No species data found in the OTU table. Skipping analysis.\n")
    quit(status = 0)  # Exit gracefully if no data
}

# Normalize the data to relative abundance
otu_data_normalized <- otu_data / sum(otu_data)

# Write species count table
write.csv(otu_data_normalized, file = species_count_file)

# Create a richness plot (example plot; you can adjust based on your needs)
richness <- rowSums(otu_data > 0)
richness_plot_data <- data.frame(Species=names(richness), Richness=richness)
png(richness_plot)
barplot(richness_plot_data$Richness, names.arg=richness_plot_data$Species, main="Richness Plot", ylab="Richness")
dev.off()

cat("Species count and richness plot generated for sample.\n")
