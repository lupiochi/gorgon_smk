rule identify_top_microbes:
    input:
        kraken2_report=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_kraken2_report.txt")
    output:
        top_nonhost=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_top_nonhost.csv")
    params:
        top_n=params["microbes"]["top_n_microbes"],
        count_threshold=params["microbes"]["count_threshold"]
    shell:
        """
        # Extract relevant species and counts, skip human species, and include only species with count above or equal to the threshold
        awk -F '\\t' '$1 ~ /s__/ {{print $1, $2}}' OFS=',' {input.kraken2_report} \
        | sed 's/.*|s__//' \
        | awk -F ',' 'NF==2 && $1 != "Homo sapiens" && $2 >= {params.count_threshold} {{print $1 "," $2}}' \
        | sort -t ',' -k2 -nr \
        | head -n {params.top_n} \
        > {output.top_nonhost}_without_header.csv

        # Check if the file contains any data before proceeding
        if [ -s {output.top_nonhost}_without_header.csv ]; then
            # Add the header to the CSV
            echo "#OTU ID,{wildcards.sample}" > {output.top_nonhost}
            cat {output.top_nonhost}_without_header.csv >> {output.top_nonhost}
        else
            # If the file is empty, create an empty CSV with just the header
            echo "#OTU ID,{wildcards.sample}" > {output.top_nonhost}
        fi

        rm -f {output.top_nonhost}_without_header.csv
        """


rule process_sample:
    input:
        csv_file=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_top_nonhost.csv")
    output:
        richness_plot=str(OUTPUT_DIR / "{sample}/figures/{sample}_{lane}_richness_plot.png"),
        species_count_file=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_species_count.csv")
    params:
        script=str(SCRIPTS_DIR / "phyloseq_analysis.r")
    conda:
        "envs/reporting_env.yaml"
    shell:
        """
        Rscript {params.script} {input.csv_file} {output.richness_plot} "NULL" {output.species_count_file}
        """

rule aggregate_samples:
    input:
        expand(str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_species_count.csv"), sample=samples["sample_id"], lane=samples["lane"])
    output:
        merged_file=str(OUTPUT_DIR / "merged_species_counts.csv")
    run:
        import pandas as pd
        
        # Read and merge all species count files
        dfs = []
        for csv_file in input:
            df = pd.read_csv(csv_file, index_col=0)
            dfs.append(df)
        
        # Concatenate all dataframes
        merged_df = pd.concat(dfs)

        # Drop any duplicate rows
        merged_df = merged_df.drop_duplicates()

        # Save the merged dataframe
        merged_df.to_csv(output.merged_file)


rule combined_analysis:
    input:
        merged_file=str(OUTPUT_DIR / "merged_species_counts.csv")
    output:
        ordination_plot=str(OUTPUT_DIR / "combined_ordination_plot.png")
    params:
        script=str(SCRIPTS_DIR / "phyloseq_analysis.r")
    conda:
        "envs/reporting_env.yaml"
    shell:
        """
        Rscript {params.script} {input.merged_file} "NULL" {output.ordination_plot} "NULL"
        """

rule multiqc:
    input:
        html_R1=expand(OUTPUT_DIR / "{sample}/fastQC/{sample}_{lane}_R1_trim_fastqc.html", sample=samples["sample_id"], lane=samples["lane"]),
        html_R2=expand(OUTPUT_DIR / "{sample}/fastQC/{sample}_{lane}_R2_trim_fastqc.html", sample=samples["sample_id"], lane=samples["lane"]),
        cutadapt_logs=expand(OUTPUT_DIR / "{sample}/{sample}_{lane}_Cutadapt.log", sample=samples["sample_id"], lane=samples["lane"])
    output:
        multiqc_report=str(OUTPUT_DIR / "{sample}/multiqc_report.html")
    conda:
        "envs/reporting_env.yaml"
    shell:
        """
        multiqc --force {OUTPUT_DIR}/{wildcards.sample} -o {OUTPUT_DIR}/{wildcards.sample}
        """

rule calculate_host_microbe_ratio:
    input:
        cram_file=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}.cram"),
        kraken2_report=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_kraken2_report.txt")
    output:
        ratio_txt=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_host_microbe_ratio.txt")
    conda:
        "envs/reporting_env.yaml"
    shell:
        """
        host_reads=$(samtools flagstat {input.cram_file} | grep 'mapped (' | head -n 1 | awk '{{print $1}}')
        microbial_reads=$(awk '{{s+=$2}} END {{print s}}' {input.kraken2_report})

        if [ "$microbial_reads" -eq 0 ]; then
            ratio="inf"
        else
            ratio=$(echo "scale=2; $host_reads / $microbial_reads" | bc)
        fi

        echo -e "Host reads: $host_reads\nMicrobial reads: $microbial_reads\nHost/Microbe Ratio: $ratio" > {output.ratio_txt}
        """



rule generate_visualizations:
    input:
        top_nonhost=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_top_nonhost.csv"),
        ratio_file=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_host_microbe_ratio.txt")
    output:
        pie_chart=str(OUTPUT_DIR / "{sample}/figures/{sample}_{lane}_host_microbe_ratio_pie_chart.png"),
        bar_chart=str(OUTPUT_DIR / "{sample}/figures/{sample}_{lane}_top_microbes_bar_chart.png")
    params:
        script=str(SCRIPTS_DIR / "visualization.py")
    conda:
        "envs/reporting_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.pie_chart})
        python {params.script} {input.top_nonhost} {input.ratio_file} {output.pie_chart} {output.bar_chart}
        """


rule compare_samples:
    input:
        species_count_files=expand(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_species_count.csv", sample=samples["sample_id"], lane=samples["lane"]),
        ratio_files=expand(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_host_microbe_ratio.txt", sample=samples["sample_id"], lane=samples["lane"])
    output:
        comparison_statistics=str(OUTPUT_DIR / "comparison_statistics.txt"),
        comparison_figures=directory(OUTPUT_DIR / "comparison_figures")
    params:
        output_dir=str(OUTPUT_DIR),
        species_count_files=lambda wildcards, input: ",".join(input.species_count_files),  # Join the input list as a comma-separated string
        ratio_files=lambda wildcards, input: ",".join(input.ratio_files)  # Join the input list as a comma-separated string
    conda:
        "envs/reporting_env.yaml"
    shell:
        """
        Rscript {SCRIPTS_DIR}/compare_samples.r {params.species_count_files} {params.ratio_files} {params.output_dir}
        """

