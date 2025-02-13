# Salmonella_Gene_Presence_Analysis_Pipeline

## Overview
This pipeline identifies the presence of a target gene across different *Salmonella enterica* serotypes using whole-genome sequencing data. It retrieves genome sequences from NCBI, processes them, and determines gene presence using alignment tools.

## Features
- **Fetches complete genome assemblies** for important *Salmonella enterica* serotypes using **NCBI Datasets** and **Entrez Direct**.
- **Processes sequencing data** with **SRA Toolkit, Trimmomatic, and SKESA** and **evaluates assemblies quality** with **SeqKit**.
- **Perform gene detection** using **BLAST+**.
- Output results in csv format, including serotype-level gene prevalence.
------
## Installation
To run this pipeline, set up a Conda environment with the required dependencies.
1. Clone the Repository
```sh
git clone https://github.com/Weifanwu66/Salmonella_Gene_Summary_Pipeline.git
cd Salmonella_Gene_Summary_Pipeline


```sh
conda env create -f environment.yml
