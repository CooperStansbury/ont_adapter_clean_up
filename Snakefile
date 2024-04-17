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
pod5_paths = os.path.abspath(config['pod5_paths'])
pod5_df = utils.load_pod5_df(pod5_paths, OUTPUT)
cell_ids = pod5_df['cell_id'].to_list() # wildcard constraints

print(f"\n======== INPUT FILES ========")
print(tabulate(pod5_df[['cell_id', 'basename']], 
      headers='keys',
      showindex=False,
      tablefmt='psql'))
  
################ ALL RULES ################
rule all:
    input:
        expand(f"{OUTPUT}pod5/{{cid}}.pod5", cid=cell_ids),
        expand(f"{OUTPUT}fast5/{{cid}}.fast5", cid=cell_ids),
        expand(f"{OUTPUT}fastq/{{cid}}.raw.fastq", cid=cell_ids),
        OUTPUT + "reports/seqkit.report.txt",
        OUTPUT + "reports/mutliqc/fastqc.html",
        OUTPUT + "reports/seqkit.porechop.report.txt",
        expand(f"{OUTPUT}quality_summary/{{cid}}.quality_summary.parquet", cid=cell_ids),
        expand(f"{OUTPUT}base_qualities/{{cid}}.base_qualities.parquet", cid=cell_ids),
        expand(f"{OUTPUT}fastqc/{{cid}}.raw_fastqc.html", cid=cell_ids),
        expand(f"{OUTPUT}porechop/{{cid}}.raw.fastq", cid=cell_ids),
        expand(f"{OUTPUT}reports/porechop_summary/{{cid}}.porechop_summary.txt", cid=cell_ids),


rule archive:
    input:
        expand(f"{OUTPUT}fastq/{{cid}}.raw.fastq.index", cid=cell_ids),
        expand(f"{OUTPUT}signal_tables/{{cid}}.raw_signal.parquet", cid=cell_ids),
        # expand(f"{OUTPUT}quality_summary_porechop/{{cid}}.quality_summary.parquet", cid=cell_ids),
        # expand(f"{OUTPUT}base_qualities_porechop/{{cid}}.base_qualities.parquet", cid=cell_ids),
        # OUTPUT + "resources/barcodes.fasta",
        # OUTPUT + "reports/mutliqc/porechop.html", 
    


rule copy_pod5:
    input:
        pod5_df['file_path'].to_list()
    output:
        protected(pod5_df['out_path'].to_list()),
    run:
        from shutil import copyfile
        for i, fpath in enumerate(input):
    
            outPath = output[i]
            copyfile(fpath, outPath)


rule get_raw_signals:
    input:
        OUTPUT + "pod5/{cid}.pod5"
    output:
        OUTPUT + "signal_tables/{cid}.raw_signal.parquet"
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """python scripts/get_raw_signal_table.py {input} {output}"""


rule make_fast5:
    input:
        OUTPUT + "pod5/{cid}.pod5"
    output:
        directory=directory(OUTPUT + "fast5/{cid}.fast5"),
        flag=touch(OUTPUT + "flags/{cid}.merge.done"),
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 4
    shell:
        """pod5 convert to_fast5 {input} \
        -t {threads} \
        --output {output.directory}"""


rule basecall:
    input:
        pod5=OUTPUT + "pod5/{cid}.pod5",
        flag=OUTPUT + "flags/{cid}.merge.done",
    output:
        protected(OUTPUT + "fastq/{cid}.raw.fastq")
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    params:
        dorado=config['dorado_path'],
        model=config['dorado_model'],
        qscore=config['dorado_min_qscore'],
    shell:
        """{params.dorado} basecaller {params.model} --emit-fastq --min-qscore {params.qscore} --no-trim {input.pod5} > {output}"""


rule f5c_index:
    input:
        fastq=OUTPUT + "fastq/{cid}.raw.fastq",
        fast5=OUTPUT + "fast5/{cid}.fast5",
    output:
        OUTPUT + "fastq/{cid}.raw.fastq.index",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    threads:
        config['threads'] // 2
    shell:
        """f5c index -t {threads} --iop 10 -d {input.fast5} {input.fastq}"""


rule get_barcode_file:
    input:
        config['barcode_path'],
    output:
        OUTPUT + "resources/barcodes.txt"
    shell:
        """cp {input} {output}"""


rule get_barcode_fasta:
    input:
        OUTPUT + "resources/barcodes.txt"
    output:
        OUTPUT + "resources/barcodes.fasta"
    shell:
        """python scripts/barcode_fasta.py {input} {output}"""


#rule locate_barcodes:
#    input:
#        fastq=OUTPUT + "fastq/{cid}.raw.fastq",
#        codes=OUTPUT + "resources/barcodes.fasta"
#    output:
#        OUTPUT + 'barcode_locations/{cid}.csv'
#    wildcard_constraints:
#        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
#    threads:
#        config['threads'] // 4
#    params:
#        m=config['n_mismatches']
#    shell:
#        """seqkit locate -m {params.m} -f {input.codes} {input.fastq} -j {threads} > {output} """


rule porechop:
    input:
        OUTPUT + "fastq/{cid}.raw.fastq",
    output:
        fastq=OUTPUT + "porechop/{cid}.raw.fastq",
        stats=OUTPUT + "porechop_stats/{cid}.porechop_stats.txt",
    threads:
        config['threads'] // 2
    shell:
        """porechop -i {input} -t {threads} -v 3 -o {output.fastq} > {output.stats}"""



rule porechop_multiqc:
    input:
        expand(f"{OUTPUT}porechop_stats/{{cid}}.porechop_stats.txt", cid=cell_ids),
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
        expand(f"{OUTPUT}fastqc/{{cid}}.raw_fastqc.html", cid=cell_ids),
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


rule quality_files:
    input:
        OUTPUT + "fastq/{cid}.raw.fastq"
    output:
        OUTPUT + "base_qualities/{cid}.base_qualities.parquet",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """python scripts/get_base_qualities.py {input} {output}"""


rule summarize_quality_files:
    input:
        OUTPUT + "base_qualities/{cid}.base_qualities.parquet",
    output:
        OUTPUT + "quality_summary/{cid}.quality_summary.parquet",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """python scripts/summarize_base_qualities.py {input} {output}"""


rule quality_files_porechop:
    input:
        OUTPUT + "porechop/{cid}.raw.fastq"
    output:
        OUTPUT + "base_qualities_porechop/{cid}.base_qualities.parquet",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """python scripts/get_base_qualities.py {input} {output}"""


rule summarize_quality_files_porechop:
    input:
        OUTPUT + "base_qualities_porechop/{cid}.base_qualities.parquet",
    output:
        OUTPUT + "quality_summary_porechop/{cid}.quality_summary.parquet",
    wildcard_constraints:
        cid='|'.join([re.escape(x) for x in set(cell_ids)]),
    shell:
        """python scripts/summarize_base_qualities.py {input} {output}"""
