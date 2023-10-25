/**
 * Rule for performing alignment of RNA-seq data using the specified alignment tool and reference genome.
 * 
 * Inputs:
 * - r1: Path to the first read file in gzipped FASTQ format.
 * - r2: Path to the second read file in gzipped FASTQ format.
 * 
 * Outputs:
 * - bam: Path to the aligned BAM file.
 * 
 * Parameters:
 * - alignment_tool_path: Path to the executable for the alignment tool to use.
 * - samtools_path: Path to the executable for samtools (optional, defaults to 'None').
 * - reference_genome_path: Path to the reference genome FASTA file.
 * 
 * Threads: Number of threads to use for alignment.
 * 
 * Resources:
 * - mem_mb: Amount of memory to allocate for the alignment process (in MB).
 * 
 * Log: Path to the log file for this rule.
 * 
 * Benchmark: Path to the benchmark file for this rule.
 * 
 * Script: Path to the Python script that performs the alignment.
 */
rule alignment:
    input:
        r1="{outdir}/01_clean_data/{sample}_val_1.fq.gz",
        r2="{outdir}/01_clean_data/{sample}_val_2.fq.gz"
    output:
        bam="{outdir}/02_aligned_data/{sample}.Aligned.sortedByCoord.out.bam"
    params:
        alignment_tool_path=config['alignment_tool_path'],
        samtools_path=config.get('samtools_path', 'None'),
        reference_genome_path=config['reference_genome_path']
    threads: 
        config["alignment_threads"]
    resources:
        mem_mb=40000
    log: "{outdir}/logs/alignment/{sample}.alignment.log"
    benchmark: "{outdir}/benchmarks/alignment/{sample}.alignment.benchmark.txt"
    script:
        "../scripts/alignment.py"