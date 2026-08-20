# Whole Exome Sequencing Pipeline

Reproducible whole-exome sequencing (WES) analysis pipeline for processing tumor and matched-normal samples, identifying somatic variants, annotating mutations, and analyzing copy-number alterations.

The workflow was developed for mouse WES data and is designed to run in an HPC environment.

## Pipeline Overview

The pipeline includes the main steps required for WES analysis:

1. Read preprocessing and quality control
2. Read mapping to the reference genome
3. Duplicate marking and base quality score recalibration
4. Somatic SNV and indel calling with GATK Mutect2
5. Variant filtering and quality assessment
6. Functional annotation with Ensembl VEP
7. Variant prioritization and visualization in R
8. Copy-number alteration analysis with CODEX2

Both Panel of Normals (PoN) and matched-normal analysis strategies are supported.

## Main Tools

* GATK / Mutect2
* Picard
* BCFtools
* Ensembl Variant Effect Predictor (VEP)
* CODEX2
* R
* Bash
* Conda
* SLURM / HPC

## Variant Analysis

Somatic variants are processed using Mutect2-derived quality metrics and VEP functional annotations.

The downstream R analysis includes:

* Variant quality filtering
* Allele-frequency analysis
* Protein-coding variant selection
* Functional consequence prioritization
* Variant summary tables
* Mutation consequence plots
* Oncoplots using `maftools`

## Copy-Number Analysis

Copy-number alterations are analyzed using CODEX2 from whole-exome sequencing coverage data.

## Configuration

Environment-specific paths, reference genomes, sample names, and computational resources are defined in the configuration file.

Example:

```bash
BASE_DIR="/path/to/WES_project/"
DATA_SHARED="/path/to/shared_data/"
CONDA_ENV_WES="/path/to/conda_envs/WESenv"
```

Users should replace these example paths with locations corresponding to their own system or HPC environment.

## Repository Structure

```text
WholeExomeSequencing-pipeline/
│
├── Scripts/
│   ├── mapping/
│   ├── variant_calling/
│   ├── annotation/
│   └── downstream_analysis/
│
├── 00_Required.sh
├── 00_conf_env.env
├── .gitignore
└── README.md
```

## Data

Raw sequencing data and large intermediate files are not included in this repository.

Files such as FASTQ, BAM, VCF, reference genomes, and analysis results should be stored separately and excluded from version control.

## Reproducibility

The workflow separates code from environment-specific configuration and uses Conda environments together with HPC resource definitions to facilitate reproducible execution across datasets and computing environments.


