"""
Rule: Calculate Host/Microbe Ratio
This rule calculates the host-to-microbe read ratio from a CRAM file and Kraken2 report.
"""
rule calculate_host_microbe_ratio:
    input:
        cram_file=str(OUTPUT_DIR / "{sample}" / "{sample}_merged.cram"),
        kraken2_report=str(OUTPUT_DIR / "{sample}" / "{sample}_kraken2_report.txt")
    output:
        ratio_txt=str(OUTPUT_DIR / "{sample}" / "{sample}_host_microbe_ratio.txt")
    conda:
        "envs/reporting_env.yaml"
    shell:
        """
        host_reads=$(samtools flagstat {input.cram_file} | grep 'mapped (' | head -n 1 | awk '{{print $1}}')
        microbial_reads=$(awk '{{s+=$2}} END {{print s}}' {input.kraken2_report})

        ratio="inf"

        if [ "$microbial_reads" -eq 0 ]; then
            ratio="inf"
        else
            ratio=$(echo "scale=2; $host_reads / $microbial_reads" | bc)
        fi

        echo -e "Host reads: $host_reads\nMicrobial reads: $microbial_reads\nHost/Microbe Ratio: $ratio" > {output.ratio_txt}
        """
