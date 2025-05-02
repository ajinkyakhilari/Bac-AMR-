# Snakefile
# ----------
configfile: "config.yaml"

import pandas as pd, os, glob

# Pull constants out of config
RAW_DIR           = config["raw_dir"]
SAMPLE_SHEET      = config["sample_sheet"]
THREADS           = config["threads"]
QUALITY           = config["quality"]
MIN_LENGTH        = config["min_length"]
ASSEMBLER         = config["assembler"].lower()
GENOME_SIZE       = config["genome_size"]
MEDAKA_MODEL      = config["medaka_model"]
MOBSUITE_DATABASE = config["mob_db_dir"]

# Map sample_id → barcode
try:
    df = pd.read_csv(SAMPLE_SHEET, dtype=str)
    if not {"barcode","sample_id"}.issubset(df.columns):
        raise KeyError
except:
    df = pd.read_csv(SAMPLE_SHEET, header=None,
                     names=["barcode","sample_id"], dtype=str)
sample2barcode = dict(zip(df.sample_id, df.barcode))
SAMPLES = list(sample2barcode)

# Contig path pattern
if ASSEMBLER == "canu":
    contig_pattern = "assembly/{sample}/canu/{sample}.contigs.fasta"
else:
    contig_pattern = "assembly/{sample}/flye/assembly.fasta"


rule all:
    input:
        expand("filtered/{sample}.fastq.gz",         sample=SAMPLES),
        expand("polished/{sample}/consensus.fasta", sample=SAMPLES),
        expand("mapping/{sample}.sorted.bam",        sample=SAMPLES),
        expand("mapping/{sample}.sorted.bam.bai",    sample=SAMPLES),
        expand("depth/{sample}.depth.txt",           sample=SAMPLES),
        expand("rgi/{sample}",                       sample=SAMPLES),
        expand("plasmids/{sample}",                  sample=SAMPLES),
        expand("quast/{sample}",                     sample=SAMPLES)


rule concat_fastq:
    input:
        lambda wc: glob.glob(os.path.join(RAW_DIR,
                                          sample2barcode[wc.sample],
                                          "*.fastq*"))
    output:
        "raw/{sample}.fastq.gz"
    conda:
        "envs/fastp.yaml"
    shell:
        "cat {input} > {output}"


rule fastp:
    input:
        "raw/{sample}.fastq.gz"
    output:
        "filtered/{sample}.fastq.gz"
    threads: THREADS
    conda:
        "envs/fastp.yaml"
    shell:
        "fastp -i {input} -o {output} --disable_adapter_trimming "
        "-w {threads} -q {QUALITY} -l {MIN_LENGTH}"


rule assembly_canu:
    input:
        "filtered/{sample}.fastq.gz"
    output:
        "assembly/{sample}/canu/{sample}.contigs.fasta"
    threads: THREADS
    conda:
        "envs/canu.yaml"
    shell:
        "canu -p {wildcards.sample} -d assembly/{wildcards.sample}/canu "
        "genomeSize={GENOME_SIZE} -nanopore {input}"


rule assembly_flye:
    input:
        "filtered/{sample}.fastq.gz"
    output:
        "assembly/{sample}/flye/assembly.fasta"
    threads: THREADS
    conda:
        "envs/flye.yaml"
    shell:
        "flye --nano-raw {input} --out-dir assembly/{wildcards.sample}/flye "
        "--threads {threads}"


rule polish:
    input:
        assembly=contig_pattern.format(sample="{sample}"),
        reads    ="filtered/{sample}.fastq.gz"
    output:
        "polished/{sample}/consensus.fasta"
    threads: THREADS
    conda:
        "envs/medaka.yaml"
    shell:
        "medaka_consensus -i {input.reads} -d {input.assembly} "
        "-m {MEDAKA_MODEL} -t {threads} "
        "-o polished/{wildcards.sample}"


rule map_reads:
    input:
        assembly="polished/{sample}/consensus.fasta",
        reads   ="filtered/{sample}.fastq.gz"
    output:
        bam="mapping/{sample}.sorted.bam",
        bai="mapping/{sample}.sorted.bam.bai"
    threads: THREADS
    conda:
        "envs/mapping.yaml"
    shell:
        """
        mkdir -p mapping
        minimap2 -ax map-ont {input.assembly} {input.reads} -t {threads} \
          | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        """


rule depth_plot:
    input:
        bam="mapping/{sample}.sorted.bam"
    output:
        depth="depth/{sample}.depth.txt"
    params:
        plot_dir="depth/plots"
    conda:
        "envs/depth.yaml"
    shell:
        """
        python scripts/depth_plot.py \
          -i {input.bam} \
          -o {output.depth} \
          -p {params.plot_dir}
        """



rule rgi:
    input:
        "polished/{sample}/consensus.fasta"
    output:
        directory("rgi/{sample}")
    conda:
        "envs/rgi.yaml"
    shell:
        """
        mkdir -p rgi/{wildcards.sample}
        rgi main -i {input} -o rgi/{wildcards.sample}/rgi_output \
          --input_type contig --clean
        """


rule mob_recon:
    input:
        fasta="polished/{sample}/consensus.fasta"
    output:
        dir=directory("plasmids/{sample}")
    log:
        "plasmids/{sample}/mob_recon.log"
    conda:
        "envs/mobsuite.yaml"
    shell:
        """
        mkdir -p {output.dir}
        mob_recon \
          -i {input.fasta} \
          -o {output.dir} \
          -d "{MOBSUITE_DATABASE}" \
          --force \
        > {log} 2>&1
        """


rule assembly_stats:
    input:
        contigs=contig_pattern.format(sample="{sample}")
    output:
        directory("quast/{sample}")
    conda:
        "envs/quast.yaml"
    shell:
        "quast -o {output} {input.contigs}"

