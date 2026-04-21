#!/bin/bash

echo "Creating Snakemake environment..."

conda create -y -n snakemake_env -c conda-forge -c bioconda snakemake

echo "Environment created."

echo "Activate environment using:"
echo "conda activate snakemake_env"

echo "Then run pipeline using:"
echo "snakemake --use-conda -j 4"
