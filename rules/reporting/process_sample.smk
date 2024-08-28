"""
Rule: Process Sample
This rule processes the top non-host microbial species CSV file and generates a richness plot and a species count file.
"""
rule process_sample:
    input:
        csv_file=str(OUTPUT_DIR / "{sample}" / "{sample}_top_nonhost.csv")
    output:
        richness_plot=str(OUTPUT_DIR / "{sample}/figures/{sample}_richness_plot.png"),
        species_count_file=str(OUTPUT_DIR / "{sample}" / "{sample}_species_count.csv")
    params:
        script=str(SCRIPTS_DIR / "process_sample_phyloseq.r")
    conda:
        "envs/reporting_env.yaml"
    shell:
        """
        num_rows=$(awk 'NR>1' {input.csv_file} | wc -l)
        
        if [ "$num_rows" -eq 0 ]; then
            echo "Skipping sample {wildcards.sample}: OTU table has no data."
            touch {output.richness_plot}
            touch {output.species_count_file}
        else
            Rscript {params.script} {input.csv_file} {output.richness_plot} {output.species_count_file}
        fi
        """
