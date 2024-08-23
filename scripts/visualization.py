import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import sys

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def ensure_directory(output_path):
    """Ensure the output directory exists."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

def plot_host_microbe_ratio_pie(file_path, output_path):
    """Generate a pie chart for host vs microbe reads."""
    try:
        logging.info(f"Reading ratio file: {file_path}")
        
        # Check if file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Ratio file not found: {file_path}")
        
        # Read and parse the ratio file
        with open(file_path, 'r') as file:
            lines = file.readlines()
            if len(lines) < 2:
                raise ValueError("Ratio file is not correctly formatted.")
            
            # Parse the host and microbial reads
            host_reads = int(lines[0].split(":")[1].strip())
            microbial_reads = int(lines[1].split(":")[1].strip())

        # Ensure the output directory exists
        ensure_directory(output_path)

        # Plot pie chart
        fig, ax = plt.subplots(figsize=(8, 8))
        labels = ['Host Reads', 'Microbial Reads']
        sizes = [host_reads, microbial_reads]
        ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140,
               colors=sns.color_palette("pastel"), pctdistance=0.85,
               wedgeprops=dict(width=0.3, edgecolor='w'))

        # Configure plot aesthetics
        ax.set_title('Host/Microbe Ratio', fontsize=16)
        ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        plt.tight_layout()

        # Save the pie chart
        plt.savefig(output_path)
        plt.close(fig)
        logging.info(f"Host/microbe ratio pie chart saved at {output_path}.")
    except Exception as e:
        logging.error(f"Error generating host/microbe ratio pie chart: {e}")
        raise

import textwrap

def plot_top_nonhost_bar(file_path, output_path):
    """Generate a vertical bar chart for the top non-host species with the sample name in the title."""
    try:
        logging.info(f"Reading top non-host file: {file_path}")
        
        # Check if file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Top non-host file not found: {file_path}")
        
        # Read the top non-host CSV file
        df_microbes = pd.read_csv(file_path)
        
        # Extract the sample name from the column header (second column name)
        sample_name = df_microbes.columns[1]
        
        # Rename the columns for easier manipulation
        df_microbes.columns = ['Species', 'Count']
        
        # Log file content for debugging
        logging.info(f"Top microbes data:\n{df_microbes.head()}")

        # Ensure the output directory exists
        ensure_directory(output_path)

        # Wrap long species names
        df_microbes['Species'] = df_microbes['Species'].apply(lambda x: '\n'.join(textwrap.wrap(x, 20)))

        # Plot vertical bar chart for the top 10 microbes
        top_nonhost = df_microbes.head(10).sort_values(by="Count", ascending=False)
        fig, ax = plt.subplots(figsize=(12, 8))  # Increased figure size for better readability
        sns.barplot(x="Species", y="Count", data=top_nonhost, palette="viridis", ax=ax)

        # Add labels on bars
        ax.bar_label(ax.containers[0], fmt='%.0f', padding=3)

        # Configure plot aesthetics
        ax.set_title(f'Top 10 Microbe Species for Sample {sample_name}', fontsize=20)
        ax.set_xlabel('Species', fontsize=16)
        ax.set_ylabel('Count', fontsize=16)
        ax.tick_params(axis='x', labelsize=14, rotation=30)  # Adjusted rotation to 30 degrees for readability
        ax.tick_params(axis='y', labelsize=14)

        # Ensure tight layout and proper aspect ratio
        plt.tight_layout()

        # Save the bar chart
        plt.savefig(output_path, dpi=300)
        plt.close(fig)
        logging.info(f"Top 10 microbe species bar chart saved at {output_path}.")
    except Exception as e:
        logging.error(f"Error generating top 10 microbe species bar chart: {e}")
        raise


def main():
    try:
        # Check that we have received exactly 4 arguments
        if len(sys.argv) != 5:
            logging.error("Incorrect number of arguments provided. Expected 4 arguments.")
            sys.exit(1)

        # Assign command line arguments to variables
        top_nonhost_file = sys.argv[1]
        ratio_file = sys.argv[2]
        pie_chart_output = sys.argv[3]
        bar_chart_output = sys.argv[4]

        logging.info(f"Received arguments: {sys.argv}")

        # Generate visualizations
        plot_host_microbe_ratio_pie(ratio_file, pie_chart_output)
        plot_top_nonhost_bar(top_nonhost_file, bar_chart_output)

    except Exception as e:
        logging.error(f"Unexpected error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
