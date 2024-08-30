library(ggplot2)
library(dplyr)
library(readr)

"""
This script generates visualizations, including:
1. A pie chart showing the ratio of host to microbe reads.
2. A bar chart displaying the top 10 non-host microbe species.

The script expects four command-line arguments:
    1. top_nonhost_file: Path to a CSV file containing the counts of non-host species.
    2. ratio_file: Path to a text file containing the ratio of host to microbe reads.
    3. pie_chart_output: Path to save the generated pie chart.
    4. bar_chart_output: Path to save the generated bar chart.
"""

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Incorrect number of arguments provided. Expected 4 arguments: <top_nonhost_file> <ratio_file> <pie_chart_output> <bar_chart_output>")
}

top_nonhost_file <- args[1]
ratio_file <- args[2]
pie_chart_output <- args[3]
bar_chart_output <- args[4]

check_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
}

ensure_directory_exists <- function(output_path) {
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
}

# Plot for host/microbe ratio pie chart
plot_host_microbe_ratio_pie <- function(file_path, output_path) {
  tryCatch({
    check_file_exists(file_path)
    
    lines <- readLines(file_path)
    if (length(lines) < 2) {
      stop("Ratio file is not correctly formatted. Something might be wrong.")
    }
    
    host_reads <- as.numeric(strsplit(lines[1], ":")[[1]][2])
    microbial_reads <- as.numeric(strsplit(lines[2], ":")[[1]][2])
    
    ensure_directory_exists(output_path)
    
    # Data for the pie chart
    pie_data <- data.frame(
      label = c("Host Reads", "Microbial Reads"),
      value = c(host_reads, microbial_reads)
    )
    
    # Plotting the pie chart
    pie_plot <- ggplot(pie_data, aes(x = "", y = value, fill = label)) +
      geom_bar(stat = "identity", width = 1, color = "white") +
      coord_polar("y") +
      theme_void() +
      theme(
        legend.title = element_blank(),
        legend.text = element_text(size = 16, color = "black"),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold", color = "#333333")
      ) +
      scale_fill_manual(values = c("#F8766D", "#00BFC4")) + # Modern contrasting colors
      geom_text(aes(label = paste0(round(value / sum(value) * 100), "%")),
                position = position_stack(vjust = 0.5), size = 6, color = "white", fontface = "bold") +
      labs(title = "Host vs Microbe Reads Distribution")
    
    ggsave(output_path, plot = pie_plot, width = 8, height = 8, dpi = 300)
    message(paste("Host/microbe ratio pie chart saved at", output_path))
    
  }, error = function(e) {
    stop(paste("Error generating host/microbe ratio pie chart:", e$message))
  })
}

# Plot for top non-host species bar chart
plot_top_nonhost_bar <- function(file_path, output_path) {
  tryCatch({
    check_file_exists(file_path)
    
    df_microbes <- read_csv(file_path)
    colnames(df_microbes) <- c("Species", "Count")
    
    ensure_directory_exists(output_path)
    
    # Select the top 10 species based on count
    top_nonhost <- df_microbes %>%
      top_n(10, Count) %>%
      arrange(desc(Count))
    
    # Plotting the bar chart with log scale for the y-axis
    bar_plot <- ggplot(top_nonhost, aes(x = reorder(Species, -Count), y = Count)) +
      geom_bar(stat = "identity", fill = "#00BFC4", color = "black", width = 0.7) +
      geom_text(aes(label = Count), vjust = -0.3, size = 4, color = "black") +
      scale_y_continuous(trans = 'log10', breaks = scales::trans_breaks("log10", function(x) 10^x)) +
      labs(title = "Top 10 Non-host Microbe Species", x = "Species", y = "Log(Count)") +
      theme_minimal(base_size = 16) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 16, face = "bold", color = "#333333"),
        axis.title.y = element_text(size = 16, face = "bold", color = "#333333"),
        plot.title = element_text(size = 20, hjust = 0.5, face = "bold", color = "#333333"),
        panel.grid.major = element_line(color = "gray80"),
        panel.grid.minor = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")
      )
    
    # Save the plot
    ggsave(output_path, plot = bar_plot, width = 12, height = 8, dpi = 300)
    message(paste("Top 10 microbe species bar chart saved at", output_path))
    
  }, error = function(e) {
    stop(paste("Error generating top 10 microbe species bar chart:", e$message))
  })
}

main <- function() {
  tryCatch({
    plot_host_microbe_ratio_pie(ratio_file, pie_chart_output)
    plot_top_nonhost_bar(top_nonhost_file, bar_chart_output)
    
  }, error = function(e) {
    stop(paste("Unexpected error occurred:", e$message))
  })
}

if (!interactive()) {
  main()
}
