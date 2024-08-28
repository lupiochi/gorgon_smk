"""
Rule: Merge CRAMs
This rule merges multiple CRAM files from different lanes into a single CRAM file, or renames the CRAM file if there is only one.
"""
rule merge_crams:
    input:
        cram_files=lambda wildcards: expand(str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}.cram"), sample=wildcards.sample, lane=SAMPLES[wildcards.sample])
    output:
        merged_cram=str(OUTPUT_DIR / "{sample}" / "{sample}_merged.cram"),
        merged_cram_index=str(OUTPUT_DIR / "{sample}" / "{sample}_merged.cram.crai")
    params:
        cores=params["compression"]["cores"]  # Use cores instead of memory for threading
    conda:
        "envs/preprocessing_env.yaml"
    shell:
        """
        num_files=$(echo {input.cram_files} | wc -w)
        
        if [ $num_files -gt 1 ]; then
            echo "Merging CRAM files..."
            samtools merge -@ {params.cores} -O CRAM {output.merged_cram} {input.cram_files}
            samtools index {output.merged_cram}
        else
            echo "Only one CRAM file, renaming..."
            mv {input.cram_files[0]} {output.merged_cram}
            mv {input.cram_files[0]}.crai {output.merged_cram_index}
        fi
        """
