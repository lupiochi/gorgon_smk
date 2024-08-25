rule index_reference_genome:
    input:
        ref_genome=str(REF_GENOME)
    output:
        ref_indexed=REF_INDEXED_FILES
    conda:
        "envs/indexing_env.yaml" 
    shell:
        """
        # Run indexing only if the index files are not found
        if [ ! -f {output[0]} ]; then
            bwa index {input.ref_genome}
            samtools faidx {input.ref_genome}  # Creates the .fai file
        fi
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

rule multiqc:
    input:
        html_R1=str(OUTPUT_DIR / "{sample}/fastQC/{sample}_{lane}_R1_trim_fastqc.html"),
        html_R2=str(OUTPUT_DIR / "{sample}/fastQC/{sample}_{lane}_R2_trim_fastqc.html"),
        cutadapt_logs=str(OUTPUT_DIR / "{sample}/{sample}_{lane}_Cutadapt.log")
    output:
        multiqc_report=str(OUTPUT_DIR / "{sample}/multiqc_report_{lane}.html")
    conda:
        "envs/reading_env.yaml"
    shell:
        """
        multiqc --force {OUTPUT_DIR}/{wildcards.sample} -o {output.multiqc_report}
        """
