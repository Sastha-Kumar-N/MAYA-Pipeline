MAYA Genome Analysis Pipeline

This pipeline performs bacterial genome analysis using Snakemake.

Tools included:
fastp
FastQC
MultiQC
SPAdes
QUAST
Prokka
KofamScan
ABRicate
MLST
Diamond
HMMER
antiSMASH
BUSCO
tRNAscan-SE
Barrnap
TRF
Jellyfish
CheckM

Installation:

1. Run installation script

bash install.sh

2. Activate environment

conda activate snakemake_env

3. Run pipeline

snakemake --use-conda -j 4
