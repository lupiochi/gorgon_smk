rule bwa_mem:
    input:
        trimmed_R1=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_R1_trim.fastq.gz"),
        trimmed_R2=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_R2_trim.fastq.gz"),
        ref_genome=str(REF_GENOME)
    output:
        aligned_bam=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}.bam")
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

rule remove_duplicates_and_compress:
    input:
        aligned_bam=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}.bam")
    output:
        cram_file=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}.cram"),
        metrics_file=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_dedup_metrics.txt"),
        cram_index=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}.cram.crai")
    params:
        memory=params["remove_duplicates"]["memory"],
        ref_genome=str(REF_GENOME)
    conda:
        "envs/preprocessing_env.yaml"
    shell:
        """
        export JAVA_OPTIONS="-Xmx{params.memory}"

        picard MarkDuplicates \
        I={input.aligned_bam} \
        O={input.aligned_bam}.dedup.bam \
        M={output.metrics_file} \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=SILENT \
        ASSUME_SORTED=true \
        QUIET=true

        samtools view -C -T {params.ref_genome} -@ 4 -o {output.cram_file} {input.aligned_bam}.dedup.bam

        samtools index {output.cram_file}

        rm {input.aligned_bam}.dedup.bam
        """

rule extract_unmapped_reads:
    input:
        cram_file=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}.cram")
    output:
        unmapped_R1=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_unmapped_1.fastq.gz"),
        unmapped_R2=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_unmapped_2.fastq.gz")
    params:
        cores=params["compression"]["cores"],
        memory=params["compression"]["memory"]
    conda:
        "envs/preprocessing_env.yaml"
    shell:
        """
        samtools view -b -f 4 {input.cram_file} | \
        samtools fastq -1 {output.unmapped_R1} -2 {output.unmapped_R2}
        
        pigz -p {params.cores} {output.unmapped_R1} {output.unmapped_R2}
        """

rule kraken2:
    input:
        unmapped_R1=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_unmapped_1.fastq.gz"),
        unmapped_R2=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_unmapped_2.fastq.gz")
    output:
        kraken2_report=str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}_kraken2_report.txt")
    params:
        kraken_db=str(KRAKEN_DB)
    conda:
        "envs/preprocessing_env.yaml"
    shell:
        """
        kraken2 --db {params.kraken_db} --paired {input.unmapped_R1} {input.unmapped_R2} \
        --report {output.kraken2_report} --use-mpa-style > /dev/null
        """
