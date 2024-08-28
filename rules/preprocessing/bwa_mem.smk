"""
Rule: BWA MEM
This rule aligns the trimmed reads to the reference genome using BWA MEM and sorts the output BAM file using Samtools.
"""
rule bwa_mem:
    input:
        trimmed_R1=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_R1_trim.fastq.gz"),
        trimmed_R2=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_R2_trim.fastq.gz"),
        ref_genome=str(REF_GENOME)
    output:
        aligned_bam=temp(str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}.bam"))
    params:
        cores=params["mapping"]["cores"],
        memory=params["mapping"]["memory"],
        rg_tag_structure=params["mapping"]["rg_tag_structure"]
    conda:
        "envs/preprocessing_env.yaml"
    shell:
        """
        bwa mem -M -t {params.cores} -R '{params.rg_tag_structure}' {input.ref_genome} {input.trimmed_R1} {input.trimmed_R2} | \
        samtools sort -@ {params.cores} -o {output.aligned_bam} 
        """
