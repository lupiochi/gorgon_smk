"""
Rule: Remove Duplicates and Compress
This rule removes duplicate reads from the BAM file and compresses the result into a CRAM file using Samtools.
"""
rule remove_duplicates_and_compress:
    input:
        aligned_bam=temp(str(OUTPUT_DIR / "{sample}" / "{sample}_{lane}.bam"))
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
