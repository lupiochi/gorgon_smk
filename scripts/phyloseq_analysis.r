# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(reshape2)

# Parse command-line arguments
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]  # This is the CSV file with species counts
ordination_plot <- args[3]  # The output ordination plot file
comparison_statistics_file <- args[4]  # Output file for comparison statistics

# Load the data
otu_data <- read.csv(data_file, row.names = 1)  # Set the first column as row names (species)

# Ensure the data has multiple samples (columns)
if (nrow(otu_data) == 0 || ncol(otu_data) < 2) {
    stop("The OTU table must have at least two samples for comparison. Please check the input data.")
}

# Calculate Pearson, Spearman correlations and R2 scores between samples
comparison_stats <- data.frame()

# Pairwise comparison between all samples
for (i in 1:(ncol(otu_data)-1)) {
  for (j in (i+1):ncol(otu_data)) {
    sample1 <- colnames(otu_data)[i]
    sample2 <- colnames(otu_data)[j]
    
    # Pearson correlation coefficient
    pcc <- cor(otu_data[, i], otu_data[, j], method = "pearson")
    
    # Spearman correlation coefficient
    scc <- cor(otu_data[, i], otu_data[, j], method = "spearman")
    
    # R2 score
    lm_model <- lm(otu_data[, i] ~ otu_data[, j])
    r2 <- summary(lm_model)$r.squared
    
    # Append the results to the dataframe
    comparison_stats <- rbind(comparison_stats, data.frame(Sample1 = sample1, Sample2 = sample2, 
                                                           Pearson_Correlation = pcc, Spearman_Correlation = scc, R2_Score = r2))
  }
}

# Write the comparison statistics to a file
write.table(comparison_stats, file = comparison_statistics_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Comparison statistics successfully written to: ", comparison_statistics_file, "\n")

# Generate heatmap of Pearson correlation coefficients between samples
correlation_matrix <- cor(otu_data, method = "pearson")

# Melt the correlation matrix for ggplot
melted_correlation_matrix <- melt(correlation_matrix)

# Generate a heatmap
p_heatmap <- ggplot(data = melted_correlation_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of Pearson Correlations Between Samples", x = "Sample", y = "Sample")

# Save the heatmap
ggsave(ordination_plot, plot = p_heatmap, width = 10, height = 6)

cat("Heatmap saved to: ", ordination_plot, "\n")
