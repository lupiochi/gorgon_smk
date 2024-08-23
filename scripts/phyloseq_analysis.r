# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(vegan)
library(data.table)

# Parse command-line arguments
args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]  # This can be a CSV for individual or merged data
richness_plot <- args[2]
ordination_plot <- args[3]
species_count_file <- args[4]  # We keep only the species count file

# Convert string "NULL" to actual NULL values in R
if (richness_plot == "NULL") richness_plot <- NULL
if (ordination_plot == "NULL") ordination_plot <- NULL
if (species_count_file == "NULL") species_count_file <- NULL

# Check if the data file is for a single sample or merged samples
is_combined_analysis <- grepl("merged", data_file)

# Load the data (either a single sample CSV or combined CSV)
otu_data <- read.csv(data_file, row.names = NULL)  # Updated: Do not set row names automatically

# Debugging: Print the dimensions of the loaded data
cat("Dimensions of the loaded OTU data: ", dim(otu_data), "\n")

# Ensure row names are unique and check data validity
if (nrow(otu_data) > 0 && ncol(otu_data) > 1) {
    # Set unique row names
    rownames(otu_data) <- make.unique(as.character(otu_data[,1]))  # Create unique row names based on the first column
    otu_data <- otu_data[,-1]  # Remove the first column since it's now used as row names
    
    # Debugging: Print the new dimensions after processing row names
    cat("Dimensions of the OTU data after processing row names: ", dim(otu_data), "\n")

    # Convert all columns to numeric, coercing non-numeric values to NA
    otu_data[] <- lapply(otu_data, function(x) as.numeric(as.character(x)))

    # Check for NA values created by coercion
    if (any(is.na(otu_data))) {
        cat("Warning: OTU data contains non-numeric values that were coerced to NA. These will be removed.\n")
        otu_data[is.na(otu_data)] <- 0  # Replace NA values with zeros or choose an appropriate strategy
    }

    # Remove rows and columns that are all zero
    otu_data <- otu_data[rowSums(otu_data) > 0, colSums(otu_data) > 0]

    # Debugging: Print the dimensions after removing zero rows and columns
    cat("Dimensions of OTU data after removing zero rows and columns: ", dim(otu_data), "\n")

    if (nrow(otu_data) == 0 || ncol(otu_data) == 0) {
        stop("OTU table has zero dimensions after filtering out zero rows and columns. Cannot proceed with analysis.")
    }

    # Convert to phyloseq OTU table
    otu_table_physeq <- otu_table(as.matrix(otu_data), taxa_are_rows = TRUE)

    # Create dummy sample metadata
    sample_metadata <- data.frame(SampleID = colnames(otu_table_physeq))
    rownames(sample_metadata) <- colnames(otu_table_physeq)
    sample_data_physeq <- sample_data(sample_metadata)

    # Combine OTU table and sample data into a phyloseq object
    physeq <- phyloseq(otu_table_physeq, sample_data_physeq)

    if (is_combined_analysis) {
        # Combined analysis with ordination
        cat("Performing combined analysis with ordination...\n")
        
        # Filter OTUs for ordination, remove OTUs present in fewer than 2 samples
        physeq <- filter_taxa(physeq, function(x) sum(x > 0) > 1, TRUE)

        # Check if there are enough samples for ordination
        if (nsamples(physeq) > 1) {
            ordination <- ordinate(physeq, method="PCoA", distance="bray")
            p2 <- plot_ordination(physeq, ordination, color="SampleID")
            if (!is.null(ordination_plot)) {
                ggsave(ordination_plot, plot=p2, width=10, height=6)
            }
        } else {
            cat("Not enough samples for ordination. Skipping.\n")
        }

    } else {
        # Single sample analysis
        cat("Performing single sample analysis...\n")

        # Calculate species counts (unique OTUs)
        species_count <- estimate_richness(physeq, split=TRUE)
        if (!is.null(species_count_file)) {
            write.csv(species_count, species_count_file, row.names=TRUE)
        }

        # Plot alpha diversity (Richness: Shannon, Simpson, Chao1)
        if (!is.null(richness_plot)) {
            p1 <- plot_richness(physeq, x="SampleID", measures=c("Shannon", "Simpson", "Chao1"))
            ggsave(richness_plot, plot=p1, width=10, height=6)
        }
    }
} else {
    # If data is invalid, stop and print a helpful error message
    stop("OTU table has zero dimensions or invalid format. Please check the input data.")
}
