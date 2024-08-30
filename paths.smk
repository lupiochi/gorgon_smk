from pathlib import Path

# Define base directories
BASE_DIR = Path("/path/to/your/project")
INPUT_DIR = Path("/path/to/your/data")
OUTPUT_DIR = BASE_DIR.joinpath("results/")  # Generic output directory for results
SCRIPTS_DIR = BASE_DIR.joinpath("scripts")  # Directory for custom scripts
KRAKEN_DB = BASE_DIR.joinpath("kraken_db/database_of_choice")  # Path to Kraken2 database
REF_GENOME_DIR = BASE_DIR.joinpath("reference_genomes")  # Directory for reference genomes

# Define reference genome and its indexed files
REF_GENOME = REF_GENOME_DIR.joinpath("reference_genome.fna")
REF_INDEXED_FILES = [
    str(REF_GENOME) + ".amb",
    str(REF_GENOME) + ".ann",
    str(REF_GENOME) + ".pac",
    str(REF_GENOME) + ".bwt",
    str(REF_GENOME) + ".fai",
    str(REF_GENOME) + ".sa"
]