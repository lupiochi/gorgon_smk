
# **GORGON: Genomics Operations for Retrieving General Organisms Notations**

<p align="center">
<img src="logo.png" alt="gorgon_logo" width="500" height="500"/>
</p>

## **Motivation**
I developed GORGON (Genomics Operations for Retrieving General Organisms Notations) as a targeted solution for researchers who need to quickly process and analyze metagenomics data (Mainly NovaSeq). While there are more comprehensive pipeline alternatives available, GORGON was specifically designed to streamline the rapid handling of large-scale datasets, enabling researchers to swiftly retrieve insights about the microorganisms present in multiple samples. The pipeline is built to facilitate the quick generation of relevant insights and visualizations, making it an efficient tool for extracting meaningful biological interpretations from complex genomic data.

GORGON automates and simplifies the intricate processes involved in metagenomic analysis, from preprocessing raw data to generating detailed reports and figures. Though it may not cover every possible analysis scenario, its focus on speed and ease of use makes it a valuable asset for timely results.

## **Key Features**
- **Modular Design**: Each step of the pipeline is modular, allowing users to customize the workflow according to their specific requirements.
- **Scalability**: GORGON is designed to handle datasets of varying sizes, from small-scale experiments to large-scale studies involving thousands of samples.
- **Reproducibility**: Built with Snakemake, the pipeline ensures reproducibility across different computational environments.
- **Comprehensive Reporting**: Generates detailed reports and visualizations to help researchers interpret their results.

## **Pipeline Overview**
### The GORGON pipeline is organized into three main stages:

#### **1. Data Reading and QC:**
- Quality control and filtering of raw sequencing data.
- Trimming and normalization of sequences.

#### **2. Data Processing:**
- Statistical analysis of genomic and microbiome data.
- Differential abundance testing.

#### **3. Reporting and Visualization:**
- Generation of summary reports.
- Visualization of results using customizable plots.

## **Environment Preparation**

### **Directory structure**
- **Remember to adjust `paths.smk` to match your work environment.**
- With `generate_samples.py`, you can easily generate the samples.tsv file to run the analysis. Make sure samples.tsv is saved in the `config` folder.
- The kraken2 database of your choice can be retrieved directly from the [kraken2 wiki](https://github.com/DerrickWood/kraken2/wiki)

```
gorgon_smk/
├── README.md
├── generate_samples.py
├── run_with_slurm.sh
├── Snakefile
├── paths.smk
├── config/
│   ├── params.yaml
│   ├── samples.tsv
├── scripts/
│   ├── compare_samples_phyloseq.r
│   ├── process_sample_phyloseq.r
│   ├── visualization.r
├── rules/
│   ├── reading/
|   |   ├── __main__.smk
|   |   ├── index_reference_genome.smk
|   |   ├── cutadapt.smk
|   |   ├── fastqc.smk
|   |   ├── multiqc.smk
│   │   ├── envs/
│   │   │   └── reading_env.yaml
│   ├── preprocessing/
|   |   ├── __main__.smk
|   |   ├── bwa_mem.smk
|   |   ├── extract_unmapped_reads.smk
|   |   ├── kraken2.smk
|   |   ├── merge_crams.smk
|   |   ├── remove_duplicates_and_compress.smk
│   │   ├── envs/
│   │   │   └── preprocessing_env.yaml
│   ├── reporting/
|   |   ├── __main__.smk
|   |   ├── aggregate_sample.smk
|   |   ├── calculate_host_microbe_ratio.smk
|   |   ├── combined_analysis.smk
|   |   ├── generate_visualizations.smk
|   |   ├── identify_top_microbes.smk
|   |   ├── process_sample.smk
│   │   ├── envs/
│   │   │   └── reporting_env.yaml
└── .gitignore
```

### **Prerequisites**
To ensure smooth execution and reproducibility of the GORGON pipeline, it is highly recommended to set up a dedicated environment. An environment, in this context, is an isolated workspace that contains all the dependencies and software required for running the pipeline without interfering with other projects or the system's global environment.

For Windows users, this can be achieved using **Windows Subsystem for Linux (WSL)**, which allows you to run a Linux distribution alongside your Windows operating system. WSL is particularly useful for running Linux-based tools and environments on a Windows machine.

### **Installing WSL on Windows**
To install WSL on your Windows machine, follow these steps:

1. Open **PowerShell** as an administrator and run:
   ```bash
   wsl --install
   ```

2. Restart your computer if prompted.

3. After the restart, set up your Linux distribution by opening the newly installed WSL terminal (e.g., Ubuntu) and following the on-screen instructions.

For more detailed instructions, refer to the official [Microsoft WSL guide](https://docs.microsoft.com/en-us/windows/wsl/install).

### **Setting Up Conda/Miniconda**
Once in a Linux-based distribution, the next step is to install **Conda/Miniconda**. Conda is an open-source package and environment management system that allows you to install and manage multiple packages and environments efficiently.

#### **Installing Miniconda**

1. Create a directory for Miniconda:
   ```bash
   mkdir -p ~/miniconda3
   ```

2. Download the Miniconda installer script:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
   ```

3. Run the installer script:
   ```bash
   bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
   ```

4. Remove the installer script to clean up:
   ```bash
   rm ~/miniconda3/miniconda.sh
   ```

5. Initialize Conda for your shell environment:
   ```bash
   ~/miniconda3/bin/conda init bash
   ~/miniconda3/bin/conda init zsh
   ```

   After running these commands, restart your terminal to apply the changes.

For further guidance on installing Miniconda, visit the [official Miniconda tutorial](https://docs.conda.io/en/latest/miniconda.html).

### **Creating a Dedicated Environment for GORGON**
After setting up Miniconda, create a dedicated Conda environment for running the GORGON pipeline. This environment will contain all the necessary dependencies and prevent conflicts with other software on your system.

1. **Create the GORGON environment:**
   ```bash
   conda create -n gorgon python=3.9
   ```

2. **Activate the environment:**
   ```bash
   conda activate gorgon
   ```

3. **Install necessary packages:**
   ```bash
   pip install snakemake pulp==2.3.1 pandas
   ```

**Note**: Mamba is often recommended as a base for Snakemake, as per [Snakemake's official tutorial](https://snakemake.readthedocs.io/en/v8.19.0/getting_started/installation.html). However, for this pipeline, running the setup above is enough.

## **Running with Slurm**
### **What is Slurm?**
Slurm (Simple Linux Utility for Resource Management) is a scalable and flexible job scheduling system commonly used in high-performance computing (HPC) environments. It manages the allocation of resources (such as CPUs, memory, and GPUs) to different jobs and queues them efficiently on the available compute nodes.

### **Running the GORGON Pipeline with Slurm**
To run the GORGON pipeline using Slurm, you may run the `run_with_slurm.sh` script with the following content, adjusting the memory and time parameters to match your settings.

```bash
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
conda activate gorgon

# Run Snakemake, using all the cores assigned by Slurm
snakemake --cores $SLURM_CPUS_PER_TASK --use-conda --latency-wait 60 --keep-going
```

### **Steps to Execute the Script**

1. Ensure that your Slurm environment is properly configured and that the necessary resources are available.

2. Submit the job script to Slurm:
   ```bash
   sbatch run_with_slurm.sh
   ```

3. Monitor the job's progress using Slurm's tools:
   ```bash
   squeue -u <your_username>
   ```

   This command will display the status of your submitted jobs.

By following these steps, you can efficiently run the GORGON pipeline in an HPC environment using Slurm, leveraging its powerful resource management capabilities.


# **Please cite me if you use this for your work**
**Luiz Felipe Piochi. GORGON: Genomics Operations for Retrieving General Organisms Notations. 2024. Retrieved from: github.com/lupiochi**