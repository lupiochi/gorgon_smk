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

        rm {input.aligned_bam}.dedup.bam
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

rule kraken2:
    input:
        unmapped_R1=str(OUTPUT_DIR / "{sample}" / "{sample}_unmapped_1.fastq.gz"),
        unmapped_R2=str(OUTPUT_DIR / "{sample}" / "{sample}_unmapped_2.fastq.gz")
    output:
        kraken2_report=str(OUTPUT_DIR / "{sample}" / "{sample}_kraken2_report.txt")
    params:
        kraken_db=str(KRAKEN_DB)
    conda:
        "envs/preprocessing_env.yaml"
    shell:
        """
        kraken2 --db {params.kraken_db} --paired {input.unmapped_R1} {input.unmapped_R2} \
        --report {output.kraken2_report} --use-mpa-style > /dev/null
        """
