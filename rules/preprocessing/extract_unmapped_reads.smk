"""
Rule: Extract Unmapped Reads
This rule extracts unmapped reads from the merged CRAM file and outputs them as FASTQ files.
"""
rule extract_unmapped_reads:
    input:
        merged_cram=str(OUTPUT_DIR / "{sample}" / "{sample}_merged.cram")
    output:
        unmapped_R1=str(OUTPUT_DIR / "{sample}" / "{sample}_unmapped_1.fastq.gz"),
        unmapped_R2=str(OUTPUT_DIR / "{sample}" / "{sample}_unmapped_2.fastq.gz")
    params:
        cores=params["compression"]["cores"],
        memory=params["compression"]["memory"]
    conda:
        "envs/preprocessing_env.yaml"
    shell:
        """
        samtools view -b -f 4 {input.merged_cram} | \
        samtools fastq -1 {output.unmapped_R1} -2 {output.unmapped_R2}
        
        pigz -p {params.cores} {output.unmapped_R1} {output.unmapped_R2}
        """
