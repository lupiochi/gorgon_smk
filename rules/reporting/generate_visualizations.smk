"""
Rule: Generate Visualizations
This rule generates a pie chart and a bar chart based on the host/microbe ratio and top microbes data.
"""
rule generate_visualizations:
    input:
        top_nonhost=str(OUTPUT_DIR / "{sample}" / "{sample}_top_nonhost.csv"),
        ratio_file=str(OUTPUT_DIR / "{sample}" / "{sample}_host_microbe_ratio.txt")
    output:
        pie_chart=str(OUTPUT_DIR / "{sample}/figures/{sample}_host_microbe_ratio_pie_chart.png"),
        bar_chart=str(OUTPUT_DIR / "{sample}/figures/{sample}_top_microbes_bar_chart.png")
    params:
        script=str(SCRIPTS_DIR / "visualization.r")
    conda:
        "envs/reporting_env.yaml"
    shell:
        """
        mkdir -p $(dirname {output.pie_chart})
        Rscript {params.script} {input.top_nonhost} {input.ratio_file} {output.pie_chart} {output.bar_chart}
        """
