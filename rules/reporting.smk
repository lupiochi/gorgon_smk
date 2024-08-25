rule identify_top_microbes:
    input:
        kraken2_report=str(OUTPUT_DIR / "{sample}" / "{sample}_kraken2_report.txt")
    output:
        top_nonhost=str(OUTPUT_DIR / "{sample}" / "{sample}_top_nonhost.csv")
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
        csv_file=str(OUTPUT_DIR / "{sample}" / "{sample}_top_nonhost.csv")
    output:
        richness_plot=str(OUTPUT_DIR / "{sample}/figures/{sample}_richness_plot.png"),
        species_count_file=str(OUTPUT_DIR / "{sample}" / "{sample}_species_count.csv")
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
        expand(str(OUTPUT_DIR / "{sample}" / "{sample}_species_count.csv"), sample=samples["sample_id"])
    output:
        merged_file=str(OUTPUT_DIR / "merged_species_counts.csv")
    run:
        import pandas as pd
        import os
        
        # Read and normalize species counts
        dfs = {}
        for csv_file in input:
            # Extract sample name from the file path (e.g., '42-i31_S23')
            sample_name = os.path.basename(csv_file).split('_species_count.csv')[0]

            # Check if sample already processed to avoid lane-specific duplicates
            if sample_name not in dfs:
                # Read each species count file
                df = pd.read_csv(csv_file, index_col=0)
                # Normalize by total counts within the sample (relative abundance)
                df = df / df.sum()
                df.columns = [sample_name]  # Rename column to sample name
                dfs[sample_name] = df  # Add dataframe to dictionary

        # Concatenate the normalized dataframes, aligning species by row index (species name)
        merged_df = pd.concat(dfs.values(), axis=1, join="outer").fillna(0)
        
        # Save the merged dataframe with each sample as a separate column
        merged_df.to_csv(output.merged_file)


rule combined_analysis:
    input:
        merged_file=str(OUTPUT_DIR / "merged_species_counts.csv")
    output:
        comparison_statistics=str(OUTPUT_DIR / "comparison_statistics_combined_analysis.txt"),
        ordination_plot=str(OUTPUT_DIR / "combined_ordination_plot.png")
    params:
        script=str(SCRIPTS_DIR / "phyloseq_analysis.r")
    conda:
        "envs/reporting_env.yaml"
    shell:
        """
        Rscript {params.script} {input.merged_file} "NULL" {output.ordination_plot} {output.comparison_statistics}
        """

rule calculate_host_microbe_ratio:
    input:
        cram_file=str(OUTPUT_DIR / "{sample}" / "{sample}_merged.cram"),
        kraken2_report=str(OUTPUT_DIR / "{sample}" / "{sample}_kraken2_report.txt")
    output:
        ratio_txt=str(OUTPUT_DIR / "{sample}" / "{sample}_host_microbe_ratio.txt")
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
        top_nonhost=str(OUTPUT_DIR / "{sample}" / "{sample}_top_nonhost.csv"),
        ratio_file=str(OUTPUT_DIR / "{sample}" / "{sample}_host_microbe_ratio.txt")
    output:
        pie_chart=str(OUTPUT_DIR / "{sample}/figures/{sample}_host_microbe_ratio_pie_chart.png"),
        bar_chart=str(OUTPUT_DIR / "{sample}/figures/{sample}_top_microbes_bar_chart.png")
    params:
        script=str(SCRIPTS_DIR / "visualization.py")
    conda:
        "envs/reporting_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.pie_chart})
        python {params.script} {input.top_nonhost} {input.ratio_file} {output.pie_chart} {output.bar_chart}
        """
