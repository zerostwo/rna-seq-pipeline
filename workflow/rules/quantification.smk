rule quantification:
    input:
        bam="{outdir}/02_aligned_data/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        counts="{outdir}/03_quantification/{sample}.counts.txt"
    log: "{outdir}/logs/quantification/{sample}.quantification.log"
    benchmark: "{outdir}/benchmarks/quantification/{sample}.quantification.benchmark.txt"
    params:
        quantification_tool_path=config['quantification_tool_path'],
        gene_annotation_path=config['gene_annotation_path']
    threads:
        config['quantification_threads']
    script:
        "../scripts/quantification.py"

rule merge_feature_counts:
    input:
        counts=expand(
            "{outdir}/03_quantification/{sample}.counts.txt", 
            sample=metadata['sample'], 
            outdir=config["outdir"]
        )
    output:
        matrix="{outdir}/03_quantification/merged_expression_matrix.txt"
    message:
        """
        Merging featureCounts results into an expression matrix.
        """
    shell:
        """
        python workflow/scripts/merge_counts.py {input.counts} > {output.matrix}
        """
