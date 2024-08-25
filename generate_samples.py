import os
import pandas as pd

def generate_samples_tsv(data_dir, output_file, extension='.fastq.gz'):
    files = os.listdir(data_dir)
    fastq_files = [f for f in files if f.endswith(extension)]
    samples_dict = {}

    for filename in fastq_files:
        parts = filename.split('_')
        if len(parts) < 5:
            continue
        
        # Join the first two parts to form the sample_id (e.g., 37_S92)
        sample_id = '_'.join(parts[:2])
        lane = parts[2]  # This will remain as the lane (e.g., L001)
        read_direction = parts[3]

        key = f"{sample_id}_{lane}"

        if key not in samples_dict:
            samples_dict[key] = {"sample_id": sample_id, "lane": lane, "forward_filename": "", "reverse_filename": ""}

        if "R1" in read_direction:
            samples_dict[key]['forward_filename'] = filename
        elif "R2" in read_direction:
            samples_dict[key]['reverse_filename'] = filename

    samples_df = pd.DataFrame.from_dict(samples_dict, orient='index').reset_index(drop=True)
    samples_df = samples_df[['sample_id', 'lane', 'forward_filename', 'reverse_filename']]
    samples_df.to_csv(output_file, sep='\t', index=False)
    print(f"Samples TSV file generated at: {output_file}")

# Your files are .fastq.gz files? If not, you can change the extension in the code below
extension = '.fastq.gz'

# Specify the directory containing FASTQ files and the output file path
data_dir = "/mnt/raw/amel_ecoapi/"
output_file = "samples.tsv"

# Generate the samples.tsv file
generate_samples_tsv(data_dir, output_file, extension)
