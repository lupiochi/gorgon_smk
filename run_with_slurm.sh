#!/bin/bash
#SBATCH --job-name=snakemake_all_samples
#SBATCH --output=snakemake_all_samples.out
#SBATCH --error=snakemake_all_samples.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G  
#SBATCH --time=96:00:00 

# Add Miniconda to PATH
export PATH=/home/usr/miniconda3/bin:$PATH

# Initialize Conda (this makes sure Conda is available in the job)
eval "$(conda shell.bash hook)"

# Activate the Conda environment that contains Snakemake
conda activate gorgon_smk

# Run Snakemake, using all the cores assigned by Slurm
snakemake --cores $SLURM_CPUS_PER_TASK --use-conda --latency-wait 60 --keep-going
