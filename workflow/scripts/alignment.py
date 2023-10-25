__author__ = "Songqi Duan"
__copyright__ = "Copyright (C) 2023 by Songqi Duan | 段松岐"
__email__ = "songqi.duan@outlook.com"
__license__ = "MIT"

from snakemake.shell import shell
import sys

alignment_tool_path = snakemake.params.alignment_tool_path
samtools_path = snakemake.params.samtools_path
input_r1 = snakemake.input.r1
input_r2 = snakemake.input.get("r2", None)
output_bam = snakemake.output.bam
reference_genome_path = snakemake.params.reference_genome_path
threads = snakemake.threads
log = snakemake.log_fmt_shell()

output_bam_prefix = output_bam.replace("Aligned.sortedByCoord.out.bam", "")

if alignment_tool_path.endswith("STAR"):
    if input_r2:
        shell(
            "{alignment_tool_path} "
            "--runThreadN {threads} "
            "--genomeDir {reference_genome_path} "
            "--readFilesIn {input_r1} {input_r2} "
            "--readFilesCommand gunzip -c "
            "--outFileNamePrefix {output_bam_prefix} "
            "--outSAMtype BAM SortedByCoordinate {log}"
        )
    else:
        shell(
            "{alignment_tool_path} "
            "--runThreadN {threads} "
            "--genomeDir {reference_genome_path} "
            "--readFilesIn {input_r1} "
            "--readFilesCommand gunzip -c "
            "--outFileNamePrefix {output_bam_prefix} "
            "--outSAMtype BAM SortedByCoordinate {log}"
        )
elif alignment_tool_path.endswith("hisat2"):
    if samtools_path is None:
        sys.exit("Error: samtools_path is required for hisat2.")
    if input_r2:
        shell(
            "{alignment_tool_path} "
            "-p {threads} "
            "-x {reference_genome_path} "
            "-1 {input_r1} "
            "-2 {input_r2} "
            "--summary-file {output_bam_prefix}summary.txt "
            "| {samtools_path} view -@ {threads} -Sbh "
            "| {samtools_path} sort -@ {threads} -o {output_bam} {log}"
        )
    else:
        shell(
            "{alignment_tool_path} "
            "-p {threads} "
            "-x {reference_genome_path} "
            "-U {input_r1} "
            "--summary-file {output_bam_prefix}summary.txt "
            "| {samtools_path} view -@ {threads} -Sbh "
            "| {samtools_path} sort -@ {threads} -o {output_bam} {log}"
        )
else:
    sys.exit("Error: Unsupported alignment tool. Currently, only 'STAR' and 'hisat2' are supported.")
