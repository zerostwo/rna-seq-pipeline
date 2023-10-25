indexed_metadata = metadata.set_index('sample')

rule quality_control:
    input:
        input_r1 = lambda wildcards: indexed_metadata.loc[wildcards.sample, 'input_r1'],
        input_r2 = lambda wildcards: indexed_metadata.loc[wildcards.sample, 'input_r2']
    output:
        output_r1 = "{outdir}/01_clean_data/{sample}_val_1.fq.gz",
        output_r2 = "{outdir}/01_clean_data/{sample}_val_2.fq.gz",
        html="{outdir}/01_clean_data/{sample}.html",
        json="{outdir}/01_clean_data/{sample}.json"
    params:
        qc_tool_path=config["qc_tool_path"]
    threads: 
        config["qc_threads"]
    log:
        "{outdir}/logs/quality_control/{sample}.quality_control.log"
    benchmark:
        "{outdir}/benchmarks/quality_control/{sample}.quality_control.benchmark.txt"
    script:
        "../scripts/quality_control.py"
