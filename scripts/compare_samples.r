# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Function to compare samples
compare_samples <- function(output_dir, species_count_files, ratio_files) {
  
  # Remove duplicate files (if any)
  species_count_files <- unique(species_count_files)
  ratio_files <- unique(ratio_files)
  
  # Aggregating data from all samples
  species_count_data <- data.frame()
  ratio_data <- data.frame()
  
  # Process each species count and ratio file
  for (i in seq_along(species_count_files)) {
    sample_name <- basename(dirname(species_count_files[i]))
    
    # Load species count data
    species_df <- read_csv(species_count_files[i])
    total_species_count <- sum(species_df[[2]])  # Assuming counts are in the second column
    
    # Load ratio data (you will need to ensure this is properly formatted and parsed)
    # This is a placeholder. Update according to the actual format of the ratio_files.
    ratio_txt <- read_lines(ratio_files[i])
    host_reads <- as.numeric(strsplit(ratio_txt[1], ":")[[1]][2])
    microbial_reads <- as.numeric(strsplit(ratio_txt[2], ":")[[1]][2])
    ratio <- host_reads / microbial_reads
    
    # Collect data for comparison
    species_count_data <- rbind(species_count_data, data.frame(sample = sample_name, total_species_count = total_species_count))
    ratio_data <- rbind(ratio_data, data.frame(sample = sample_name, host_reads = host_reads, microbial_reads = microbial_reads, host_microbe_ratio = ratio))
  }
  
  # Inspect the data before plotting
  print(species_count_data)
  
  # Prepare directories for output
  comparison_dir <- file.path(output_dir, "comparison_figures")
  dir.create(comparison_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Perform statistical tests between pairs of samples (t-test)
  comparison_results <- data.frame()
  samples <- unique(species_count_data$sample)
  
  for (i in seq_along(samples)) {
    for (j in seq_along(samples)) {
      if (i < j) {
        sample1_data <- species_count_data %>% filter(sample == samples[i])
        sample2_data <- species_count_data %>% filter(sample == samples[j])
        
        # Check if variance exists and is not NA in both groups before running the t-test
        sample1_var <- var(sample1_data$total_species_count, na.rm = TRUE)
        sample2_var <- var(sample2_data$total_species_count, na.rm = TRUE)
        
        if (!is.na(sample1_var) && !is.na(sample2_var) && sample1_var > 0 && sample2_var > 0) {
          # Perform t-test if there is variation
          test_result <- t.test(sample1_data$total_species_count, sample2_data$total_species_count)
          comparison_results <- rbind(comparison_results, data.frame(
            sample1 = samples[i],
            sample2 = samples[j],
            statistic = test_result$statistic,
            p_value = test_result$p.value
          ))
        } else {
          # Log that the comparison could not be performed due to lack of variance or NA variance
          comparison_results <- rbind(comparison_results, data.frame(
            sample1 = samples[i],
            sample2 = samples[j],
            statistic = NA,
            p_value = NA
          ))
          cat("No valid variance between", samples[i], "and", samples[j], "skipping t-test.\n")
        }
      }
    }
  }
  
  # Write comparison results to file
  write_csv(comparison_results, file.path(output_dir, "comparison_statistics.txt"))
  
  # Create visualizations with solid background (no transparency)
  # 1. Total species count comparison
  p1 <- ggplot(species_count_data, aes(x = sample, y = total_species_count)) +
    geom_boxplot(outlier.colour = "red", outlier.shape = 8) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = "black")) +
    ggtitle("Total Species Count Across Samples")
  
  ggsave(file.path(comparison_dir, "total_species_count_comparison.png"), plot = p1, width = 10, height = 6, bg = "white")
  
  # 2. Host/Microbe ratio comparison
  p2 <- ggplot(ratio_data, aes(x = sample, y = host_microbe_ratio)) +
    geom_boxplot(outlier.colour = "red", outlier.shape = 8) +
    theme_minimal() +
    theme(panel.background = element_rect(fill = "white", color = "black")) +
    ggtitle("Host/Microbe Ratio Across Samples")
  
  ggsave(file.path(comparison_dir, "host_microbe_ratio_comparison.png"), plot = p2, width = 10, height = 6, bg = "white")
  
  cat("Comparative analysis and visualizations generated successfully.\n")
}

# Main script execution
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Expected at least 3 arguments: species_count_files, ratio_files, and output_dir")
}

species_count_files <- unlist(strsplit(args[1], ","))
ratio_files <- unlist(strsplit(args[2], ","))
output_dir <- args[3]

compare_samples(output_dir, species_count_files, ratio_files)
