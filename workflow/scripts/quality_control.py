__author__ = "Songqi Duan"
__copyright__ = "Copyright (C) 2023 by Songqi Duan | 段松岐"
__email__ = "songqi.duan@outlook.com"
__license__ = "MIT"

from snakemake.shell import shell
import sys

qc_tool_path = snakemake.params.qc_tool_path
input_r1 = snakemake.input.input_r1
input_r2 = snakemake.input.get("input_r2", None)
output_r1 = snakemake.output.output_r1
output_r2 = snakemake.output.get("output_r2", None)
threads = snakemake.threads
json = snakemake.output.get("json", "")
html = snakemake.output.get("html", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if qc_tool_path.endswith("fastp"):
    if input_r2:
        shell(
            "{qc_tool_path} --in1 {input_r1} --in2 {input_r2} "
            "--out1 {output_r1} --out2 {output_r2} "
            "--thread {threads} "
            "--json {json} "
            "--html {html} {log}"
        )
    else:
        shell(
            "{qc_tool_path} --in1 {input_r1} "
            "--out1 {output_r1} "
            "--thread {threads} "
            "--json {json} "
            "--html {html} {log}"
        )
else:
    sys.exit("Error: Unsupported quality control tool. Currently, only 'fastp' is supported.")
