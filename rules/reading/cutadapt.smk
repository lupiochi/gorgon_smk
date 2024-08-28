"""
Rule: Cutadapt
This rule trims adapters and filters low-quality bases from sequencing reads using Cutadapt.
"""
rule cutadapt:
    input:
        R1=str(INPUT_DIR / "{sample}_{lane}_R1_001.fastq.gz"),
        R2=str(INPUT_DIR / "{sample}_{lane}_R2_001.fastq.gz")
    output:
        trimmed_R1=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_R1_trim.fastq.gz"),
        trimmed_R2=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_R2_trim.fastq.gz"),
        log_file=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_Cutadapt.log")
    params:
        cores=params["cutadapt"]["cores"],
        memory=params["cutadapt"]["memory"],
        adapter1=params["cutadapt"]["adapter1"],
        adapter2=params["cutadapt"]["adapter2"]
    conda:
        "envs/reading_env.yaml"
    shell:
        """
        cutadapt --cores {params.cores} -q 20 -m 20 -a {params.adapter1} -A {params.adapter2} \
        --max-n 0 -o {output.trimmed_R1} -p {output.trimmed_R2} {input.R1} {input.R2} > {output.log_file}
        """
