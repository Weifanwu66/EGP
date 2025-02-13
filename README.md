# Salmonella_Gene_Summary_Pipeline

## Overview
This pipeline identifies the presence of a target gene across different *Salmonella enterica* serotypes using whole-genome sequencing data. It retrieves genome sequences from NCBI, processes them, and determines gene presence using alignment tools.

## Features
- **Fetches complete genome assemblies** for important *Salmonella enterica* serotypes using **NCBI Datasets** and **Entrez Direct**.
- **Fetches and processes raw sequencing data** with **SRA Toolkit, Trimmomatic, and SKESA** and **evaluates assemblies quality** with **SeqKit**.
- **Perform gene detection** using **BLAST+**.
- Supports both **Heavy-weight and Light-weight Modes**:
  - Light-weight Mode: Uses only pre-exisiting complete genomes retrieved from NCBI database for analysis.
  - Heavy-weight Mode: While still analyzing complete genomes, if a serotype has fewer than 30 complete genomes, the pipeline will automatically retrieve SRA data, assemble the raw sequences, and assess the quality of assemblies.
- Output results in csv format, including serotype-level gene prevalence.
------
## Installation
To run this pipeline, set up a Conda environment with the required dependencies.
1. Clone the Repository
```sh
git clone https://github.com/Weifanwu66/Salmonella_Gene_Summary_Pipeline.git
cd Salmonella_Gene_Summary_Pipeline
```
2. Create and Activate the Conda Environment
The pipeline requires a Conda environment with all necessary dependencies. To create and activate it, run:
```sh
conda env create -f environment.yml
conda activate Salmonella_Gene_Summary
```
To verify the installation, check if all tools are installed:
```sh
conda list | grep -E "sra-tools|seqkit|trimmomatic|skesa|blast|ncbi-datasets-cli|entrez-direct"
```
If any package is missing, please install it manually:
```sh
conda install -c bioconda <package_name>
```
-----
## Dependencies
1. [NCBI Datasets] https://github.com/ncbi/datasets
-----
## Flags and Options
The pipeline provides several flgas to customize execution:
- -g : 
