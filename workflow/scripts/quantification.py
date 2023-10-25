__author__ = "Songqi Duan"
__copyright__ = "Copyright (C) 2023 by Songqi Duan | 段松岐"
__email__ = "songqi.duan@outlook.com"
__license__ = "MIT"

from snakemake.shell import shell
import sys

quantification_tool_path = snakemake.params.quantification_tool_path
bam_file = snakemake.input.bam
output_counts = snakemake.output.counts
gene_annotation_path = snakemake.params.gene_annotation_path
threads = snakemake.threads
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if quantification_tool_path.endswith("featureCounts"):
    shell(
        "{quantification_tool_path} "
        "-a {gene_annotation_path} "
        "-g gene_name "
        "-o {output_counts} "
        "-T {threads} "
        "-p {bam_file} {log}"
    )
elif quantification_tool_path.endswith("htseq-count"):
    shell(
        "{quantification_tool_path} "
        "-f bam "
        "-r pos "
        "-t gene "
        "-i gene_id {bam_file} {gene_annotation_path} > {output_counts} {log}"
    )
else:
    sys.exit("Error: Unsupported quantification tool. Only 'featureCounts' and 'htseq-count' are supported.")
