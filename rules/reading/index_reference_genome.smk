"""
Rule: Index Reference Genome
This rule indexes the reference genome using BWA and Samtools.
"""
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
