import yaml
import pandas as pd
from pathlib import Path

# Load configurations
config_path = Path("config")
params = yaml.safe_load(open(config_path / "params.yaml"))
samples = pd.read_table(config_path / "samples.tsv", comment="#", dtype="str")

# Define the SAMPLES variable from the samples DataFrame
SAMPLES = samples.groupby("sample_id")["lane"].apply(list).to_dict()

# Sanity check for samples file
assert "sample_id" in samples.columns, "Missing 'sample_id' column in samples.tsv"

# Include paths and rule files
include: "paths.smk"
include: "rules/reading.smk"
include: "rules/preprocessing.smk"
include: "rules/reporting.smk"

# Final 'all' rule to indicate what needs to be completed
rule all:
    input:
        # These are individual outputs that will be generated for each sample
        expand(OUTPUT_DIR / "{sample}/figures/{sample}_host_microbe_ratio_pie_chart.png", sample=samples["sample_id"]),
        expand(OUTPUT_DIR / "{sample}/figures/{sample}_top_microbes_bar_chart.png", sample=samples["sample_id"]),
        expand(OUTPUT_DIR / "{sample}/figures/{sample}_richness_plot.png", sample=samples["sample_id"]),
        expand(OUTPUT_DIR / "{sample}" / "multiqc_report_{lane}.html", sample=samples["sample_id"], lane=samples["lane"]), 
        str(OUTPUT_DIR / "combined_ordination_plot.png")