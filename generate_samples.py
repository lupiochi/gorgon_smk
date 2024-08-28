import os
import pandas as pd

# User-configurable settings and examples
"""
# Example 1 (Common Illumina setup with default R1/R2 and underscore separation):
data_dir = "/path/to/your/data/"
output_file = "samples.tsv"
extension = '.fastq.gz'
separator = '_'
sample_id_split_indices = [0, 1]  # Parts to join for sample_id (e.g., 42-i3_S9)
lane_index = 2  # Position of the lane in the filename (e.g., L001)
read_direction_indicator = "R1/R2"  # Indicator for read direction in the filename

# Example 2 (Alternative Illumina setup with hyphen separation):
data_dir = "/path/to/your/data/"
output_file = "samples.tsv"
extension = '.fq.gz'
separator = '-'
sample_id_split_indices = [0, 1]  # Parts to join for sample_id (e.g., 42-i3-S9)
lane_index = 3  # Position of the lane in the filename (e.g., L001)
read_direction_indicator = "R1/R2"  # Indicator for read direction in the filename
"""

# User settings for customization
data_dir = "/mnt/raw/amel_ecoapi/"
output_file = "samples.tsv"
extension = '.fastq.gz'
separator = '_'
sample_id_split_indices = [0, 1]
lane_index = 2
read_direction_indicator = "R1/R2"

def generate_samples_tsv(data_dir, output_file, extension, separator, sample_id_split_indices, lane_index, read_direction_indicator):
    files = os.listdir(data_dir)
    fastq_files = [f for f in files if f.endswith(extension)]
    samples_dict = {}

    for filename in fastq_files:
        parts = filename.split(separator)
        if len(parts) < max(sample_id_split_indices) + 1:
            continue
        
        sample_id = separator.join([parts[i] for i in sample_id_split_indices])
        lane = parts[lane_index]
        read_direction = parts[lane_index + 1]

        key = f"{sample_id}_{lane}"

        if key not in samples_dict:
            samples_dict[key] = {"sample_id": sample_id, "lane": lane, "forward_filename": "", "reverse_filename": ""}

        if read_direction_indicator.split('/')[0] in read_direction:
            samples_dict[key]['forward_filename'] = filename
        elif read_direction_indicator.split('/')[1] in read_direction:
            samples_dict[key]['reverse_filename'] = filename

    samples_df = pd.DataFrame.from_dict(samples_dict, orient='index').reset_index(drop=True)
    samples_df = samples_df[['sample_id', 'lane', 'forward_filename', 'reverse_filename']]
    samples_df.to_csv(output_file, sep='\t', index=False)
    print(f"Samples TSV file generated at: {output_file}")

generate_samples_tsv(data_dir, output_file, extension, separator, sample_id_split_indices, lane_index, read_direction_indicator)
