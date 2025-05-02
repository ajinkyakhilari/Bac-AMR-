# Bac-AMR
The Bacterial Whole Genome Analysis Pipeline 
This repository contains a Snakemake-based pipeline for assembly, polishing, annotation, and plasmid reconstruction of nanopore-sequenced bacterial isolates.

### Pipeline Overview

At a high level, the Bac‑AMR pipeline automates the end‑to‑end processing of nanopore‑sequenced bacterial isolates. Starting from raw per‑barcode FASTQ files, it performs quality filtering (fastp), de novo assembly (Flye or Canu), consensus polishing (Medaka), read mapping and coverage calculation (minimap2 + samtools), antimicrobial resistance gene identification (RGI with the CARD database), plasmid reconstruction (MOB‑suite), and final assembly metrics reporting (QUAST). Each step is encapsulated in a Snakemake rule with its own Conda environment, ensuring reproducibility and scalable parallel execution.

## Directory Structure

```bash
├── Snakefile
├── config.yaml
├── envs/
│   ├── fastp.yaml
│   ├── canu.yaml
│   ├── flye.yaml
│   ├── medaka.yaml
│   ├── mapping.yaml
│   ├── depth.yaml
│   ├── rgi.yaml
│   ├── mobsuite.yaml
│   └── quast.yaml
├── scripts/
│   └── depth_plot.py
└── README.md
```

#### Snakefile: Workflow definition.

#### config.yaml: User-configurable parameters (paths, threads, quality thresholds, etc.).

#### envs/: Conda environment specifications for each rule.

#### scripts/: Auxiliary scripts (e.g., depth_plot.py for coverage plotting).

## Installation

#### Clone the repo

```bash
git clone https://github.com/your-org/Bac-AMR.git
cd Bac-AMR
```

#### Install Snakemake (if not already available)

Using Conda:

```bash
conda create -n snakemake -c conda-forge snakemake
conda activate snakemake
```

Or via pip:

```bash
pip install snakemake
```

#### Install RGI and MOB-suite

```bash
conda install -c bioconda -c conda-forge rgi mobsuite
```

#### Download CARD database and variant detection models (for RGI)

```bash
rgi load \
  --card_json https://card.mcmaster.ca/latest/data \
  --model_json https://card.mcmaster.ca/latest/models \
  --local
```

#### Download MOB-suite databases

```bash
mob_init --db-dir /path/to/mob_suite/databases
```

#### Prepare environments

The first run of Snakemake will automatically create the Conda environments under envs/. Ensure --use-conda is enabled.

#### Configuration

Edit config.yaml to set project paths and parameters. Example:

```bash
raw_dir: "/data/fastq"
sample_sheet: "samples.csv"
threads: 16
quality: 10
min_length: 1000
assembler: "flye"
genome_size: "5m"
medaka_model: "r941_prom_sup_g507"
mob_db_dir: "/path/to/mob_suite/databases"
```

- raw_dir: directory containing per-barcode FASTQ subfolders.

- sample_sheet: CSV mapping barcode,sample_id.

- threads: number of threads to use.

- quality, min_length: fastp filtering thresholds.

- assembler: flye or canu.

- genome_size: genome size parameter for Canu.

- medaka_model: model for Medaka polishing.

- mob_db_dir: path to MOB-suite database directory.

#### Usage

Basic run

From within the project directory:

```bash
snakemake --cores 16 --use-conda
```

This will consume all samples defined in config.yaml and produce outputs in ```raw/```, ```filtered/```, ```assembly/```, ```polished/```, ```mapping/```, ```depth/```, ```rgi/```,``` plasmids/```, and ```quast/```

### Customizing on the fly

You can override config parameters directly from the command line using ```--config``` For example:

```bash
snakemake \
  --cores 8 \
  --use-conda \
  --config raw_dir=/new/raw/path \
    threads=8 \
    assembler=canu
```

#### Dry-run and DAG

Preview what will be executed without running:

```bash
snakemake -n --cores 16 --use-conda
```

#### Plot the workflow DAG:

```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

### Running from a different directory

If you’re outside the project folder, specify both the Snakefile and config file:

```bash
snakemake \
  --snakefile /path/to/Bac-AMR/Snakefile \
  --configfile /path/to/Bac-AMR/config.yaml \
  --cores 4 \
  --use-conda
  --schedular greedy
```

### Depth Plotting

```Coverage depth is computed and plotted by scripts/depth_plot.py. The PNGs are saved in depth/plots.```

Contact

For issues or questions, please open an issue on GitHub or contact the maintainers.



