from pathlib import Path

# Define base directories
BASE_DIR = Path("/home/lpiochi/")
INPUT_DIR = Path("/mnt/raw/amel_ecoapi")
OUTPUT_DIR = Path("/mnt/share/bees/ecoapi/results_piochi_smk")
SCRIPTS_DIR = BASE_DIR.joinpath("snake_v4/scripts/")
KRAKEN_DB = BASE_DIR.joinpath("minikraken2_db/minikraken2_v2_8GB_201904_UPDATE")
REF_GENOME_DIR = Path("/mnt/share/bees/RefGenomes/Amel3.1/")

# Define reference genome and its indexed files
REF_GENOME = REF_GENOME_DIR.joinpath("Amel3.1.fna")
REF_INDEXED_FILES = [
    str(REF_GENOME) + ".amb",
    str(REF_GENOME) + ".ann",
    str(REF_GENOME) + ".pac",
    str(REF_GENOME) + ".bwt",
    str(REF_GENOME) + ".fai",
    str(REF_GENOME) + ".sa"
]