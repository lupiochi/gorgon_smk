"""
Rule: Identify Top Microbes
This rule extracts the top non-host microbial species from the Kraken2 report, filtering by a count threshold.
"""
rule identify_top_microbes:
    input:
        kraken2_report=str(OUTPUT_DIR / "{sample}" / "{sample}_kraken2_report.txt")
    output:
        top_nonhost=str(OUTPUT_DIR / "{sample}" / "{sample}_top_nonhost.csv")
    params:
        top_n=params["microbes"]["top_n_microbes"],
        count_threshold=params["microbes"]["count_threshold"]
    shell:
        """
        awk -F '\\t' '$1 ~ /s__/ {{print $1, $2}}' OFS=',' {input.kraken2_report} \
        | sed 's/.*|s__//' \
        | awk -F ',' 'NF==2 && $1 != "Homo sapiens" && $2 >= {params.count_threshold} {{print $1 "," $2}}' \
        | sort -t ',' -k2 -nr \
        | head -n {params.top_n} \
        > {output.top_nonhost}_without_header.csv

        if [ -s {output.top_nonhost}_without_header.csv ]; then
            echo "#OTU ID,{wildcards.sample}" > {output.top_nonhost}
            cat {output.top_nonhost}_without_header.csv >> {output.top_nonhost}
        else
            echo "#OTU ID,{wildcards.sample}" > {output.top_nonhost}
        fi

        rm -f {output.top_nonhost}_without_header.csv
        """
