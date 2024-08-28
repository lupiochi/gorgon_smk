library(phyloseq)
library(ggplot2)
library(reshape2)
library(ggrepel)  # To handle sample label overlap

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop("Insufficient arguments provided. Usage: Rscript script.R <data_file> <ordination_plot> <comparison_statistics_file>")
}
data_file <- args[1]
ordination_plot <- args[2]
comparison_statistics_file <- args[3]

# Load the data and check format
otu_data <- read.csv(data_file, row.names = 1)
# Check if there are at least two samples (columns) for comparison
if (ncol(otu_data) < 2) {
    cat("The OTU table must have at least two samples for comparison. Skipping analysis.\n")
    quit(status = 0)  # Exit gracefully
}

# Continue with analysis if there are two or more samples
comparison_stats <- data.frame()

# Pairwise comparison between samples
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

write.table(comparison_stats, file = comparison_statistics_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Comparison statistics successfully written to: ", comparison_statistics_file, "\n")

# Ordination Analysis (PCA)
otu_data_t <- t(otu_data)  # Transpose OTU data to have samples as rows and taxa as columns
otu_pca <- prcomp(otu_data_t, scale. = TRUE)

# Create a data frame for plotting
ordination_df <- as.data.frame(otu_pca$x)
ordination_df$Sample <- rownames(ordination_df)

# Set theme for the plot
theme_publication <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_line(color = "gray90"),
  panel.grid.minor = element_blank(),
  axis.line = element_line(color = "black"),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12),
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
)

# Define a threshold for labeling
max_samples_to_label <- 50  # Set a threshold for maximum samples to label

# Plotting with sample labels or without, based on the number of samples
if (nrow(ordination_df) <= max_samples_to_label) {
  ordination_plot_obj <- ggplot(ordination_df, aes(x = PC1, y = PC2, label = Sample)) +
    geom_point(size = 4, color = "blue") +
    geom_text_repel(size = 3.5) +  # Using ggrepel to avoid label overlap
    labs(title = "PCA Ordination Plot", x = "PC1", y = "PC2") +
    theme_publication
} else {
  ordination_plot_obj <- ggplot(ordination_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Sample), size = 4) +
    labs(title = "PCA Ordination Plot", x = "PC1", y = "PC2") +
    theme_publication +
    theme(legend.position = "right")  # Display a legend instead of labels
}

# Save the ordination plot
ggsave(ordination_plot, plot = ordination_plot_obj, width = 10, height = 7)
cat("Ordination plot successfully saved at: ", ordination_plot, "\n")
