"""
Rule: Aggregate Samples
This rule aggregates the species count data from multiple samples and generates a merged species count CSV file.
"""
rule aggregate_samples:
    input:
        expand(str(OUTPUT_DIR / "{sample}" / "{sample}_species_count.csv"), sample=samples["sample_id"])
    output:
        merged_file=str(OUTPUT_DIR / "merged_species_counts.csv")
    run:
        import pandas as pd
        import os
        
        dfs = {}
        for csv_file in input:
            sample_name = os.path.basename(csv_file).split('_species_count.csv')[0]

            if sample_name not in dfs:
                df = pd.read_csv(csv_file, index_col=0)
                df = df / df.sum()
                df.columns = [sample_name]
                dfs[sample_name] = df

        merged_df = pd.concat(dfs.values(), axis=1, join="outer").fillna(0)
        merged_df.to_csv(output.merged_file)
