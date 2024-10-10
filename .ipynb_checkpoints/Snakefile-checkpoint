import pandas as pd
import yaml
from pathlib import Path
import re
import os
import sys
from utils import utils
from snakemake.utils import Paramspace
from tabulate import tabulate

BASE_DIR = Path(workflow.basedir)
configfile: str(BASE_DIR) + "/config/config.yaml"
OUTPUT = config['output_path']

# load the basecalls
fastq_paths = os.path.abspath(config['fastq_paths'])
fastq_df = utils.load_fastq_df(fastq_paths, OUTPUT)
cell_ids = fastq_df['cell_id'].to_list() # wildcard constraints

print(f"\n======== INPUT FILES ========")
print(tabulate(fastq_df[['cell_id', 'basename']], 
      headers='keys',
      showindex=False,
      tablefmt='psql'))
  
################ ALL RULES ################
rule all:
    input:
        expand(OUTPUT + "fastq/{cid}.raw.fastq", cid=cell_ids),
        OUTPUT + "reports/seqkit.report.txt",
        OUTPUT + "reports/seqkit.porechop.report.txt",
        expand(OUTPUT  + "porechop/{cid}.raw.fastq", cid=cell_ids),
        expand(OUTPUT  + "reports/porechop_summary/{cid}.porechop_summary.txt", cid=cell_ids),
        OUTPUT + "reports/adapter_summary.csv",
        
        
rule fastqc:
    input:
        expand(OUTPUT  + "fastqc/{cid}.raw_fastqc.html", cid=cell_ids),
        OUTPUT + "reports/mutliqc/fastqc.html",


rule copy_fastq:
    input:
        fastq_df['file_path'].to_list()
    output:
        fastq_df['out_path'].to_list(),
    run:
        from shutil import copyfile
        for i, fpath in enumerate(input):
    
            outPath = output[i]
            copyfile(fpath, outPath)


rule get_barcode_file:
    input:
        config['barcode_path'],
    output:
        OUTPUT + "resources/barcodes.txt"
    shell:
        """cp {input} {output}"""


rule porechop:
    input:
        OUTPUT + "fastq/{cid}.raw.fastq",
    output:
        fastq=OUTPUT + "porechop/{cid}.raw.fastq",
        stats=OUTPUT + "porechop_stats/{cid}.porechop_stats.txt",
    threads:
        config['threads'] // 2
    conda:
        "porechop"
    shell:
        """porechop -i {input} -t {threads} -v 3 -o {output.fastq} > {output.stats}"""


rule porechop_multiqc:
    input:
        expand(OUTPUT + "porechop_stats/{cid}.porechop_stats.txt", cid=cell_ids),
    output:
        OUTPUT + "reports/mutliqc/porechop.html",
        directory(OUTPUT + "reports/mutliqc/porechop_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.5.2/bio/multiqc"


rule fastqc_report:
    input:
        OUTPUT + "fastq/{cid}.raw.fastq",
    output:
        html=OUTPUT + "fastqc/{cid}.raw_fastqc.html",
        zip=OUTPUT + "fastqc/{cid}.raw_fastqc.zip" 
    params: f"--quiet --contaminants {config['fastqc_contaminants']} --adapters {config['fastqc_adapters']}"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads: 
        config['threads'] // 4
    wrapper:
        "v1.14.1/bio/fastqc"


rule fastqc_multiqc:
    input:
        expand(OUTPUT + "fastqc/{cid}.raw_fastqc.html", cid=cell_ids),
    output:
        OUTPUT + "reports/mutliqc/fastqc.html",
        directory(OUTPUT + "reports/mutliqc/fastqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.5.2/bio/multiqc"


rule fastq_report:
    input:
        expand(OUTPUT + "fastq/{cid}.raw.fastq", cid=cell_ids),
    output:
        OUTPUT + "reports/seqkit.report.txt",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4
    conda:
        "porechop"
    shell:
        """seqkit stats -a -b -j {threads} {input} -o {output}"""


rule porechop_report:
    input:
        expand(OUTPUT + "porechop/{cid}.raw.fastq", cid=cell_ids),
    output:
        OUTPUT + "reports/seqkit.porechop.report.txt",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4
    conda:
        "porechop"
    shell:
        """seqkit stats -a -b -j {threads} {input} -o {output}"""


rule porechop_summary:
    input:
        OUTPUT + "porechop_stats/{cid}.porechop_stats.txt",
    output:
        OUTPUT + "reports/porechop_summary/{cid}.porechop_summary.txt",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """cat {input} | grep 'reads' | grep 'adapters' | grep -v 'containing' > {output}"""


rule adapter_content_summary:
    input:
        OUTPUT + "porechop_stats/{cid}.porechop_stats.txt",
    output:
        OUTPUT + "reports/adapter_content/{cid}.txt"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """cat {input} | grep 'read coords' > {output}"""



rule adapter_summary:
    input:
        expand(OUTPUT + "reports/adapter_content/{cid}.txt", cid=cell_ids),
    output:
        OUTPUT + "reports/adapter_summary.csv"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    conda:
        "bioinf"
    shell:
        """python scripts/compile_adapter_summary.py {output} {input}"""
