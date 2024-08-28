"""
Rule: FastQC
This rule performs quality control checks on the trimmed reads using FastQC.
"""
rule fastqc:
    input:
        trimmed_R1=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_R1_trim.fastq.gz"),
        trimmed_R2=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_R2_trim.fastq.gz")
    output:
        html_R1=str(OUTPUT_DIR / "{sample}/fastQC/{sample}_{lane}_R1_trim_fastqc.html"),
        html_R2=str(OUTPUT_DIR / "{sample}/fastQC/{sample}_{lane}_R2_trim_fastqc.html")
    params:
        cores=params["fastqc"]["cores"],
        memory=params["fastqc"]["memory"]
    conda:
        "envs/reading_env.yaml"
    shell:
        """
        fastqc --threads {params.cores} -o {OUTPUT_DIR}/{wildcards.sample}/fastQC {input.trimmed_R1} {input.trimmed_R2}
        """
