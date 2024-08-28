"""
Rule: Combined Analysis
This rule performs combined analysis on all samples and generates comparison statistics and an ordination plot.
"""
rule combined_analysis:
    input:
        merged_file=str(OUTPUT_DIR / "merged_species_counts.csv")
    output:
        ordination_plot=str(OUTPUT_DIR / "combined_ordination_plot.png"),
        comparison_statistics_file=str(OUTPUT_DIR / "comparison_statistics.tsv")
    params:
        script=str(SCRIPTS_DIR / "compare_samples_phyloseq.r")
    conda:
        "envs/reporting_env.yaml"
    shell:
        """
        Rscript {params.script} {input.merged_file} {output.ordination_plot} {output.comparison_statistics_file}
        """
