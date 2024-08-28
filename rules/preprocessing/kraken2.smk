"""
Rule: Kraken2
This rule runs Kraken2 to classify the unmapped reads against a Kraken2 database.
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
