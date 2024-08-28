"""
Rule: MultiQC
This rule aggregates results from Cutadapt and FastQC into a summary report using MultiQC.
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
