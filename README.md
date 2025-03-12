# Enterobacteriaceae_Gene_Summary_Pipeline

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
git clone https://github.com/Weifanwu66/SGP.git
cd SGP
```
2. Create and Activate the Conda Environment
The pipeline requires a Conda environment with all necessary dependencies. To create and activate it, run:
```sh
conda env create -f environment.yml
conda activate SGP
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
1. NCBI Datasets: https://github.com/ncbi/datasets
2. Entrez Direct: https://www.ncbi.nlm.nih.gov/books/NBK179288/
3. NCBI BLAST: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html
4. SRA Toolkit: https://github.com/ncbi/sra-tools/wiki
5. Trimmomatic: https://github.com/usadellab/Trimmomatic/releases
6. SKESA: https://github.com/ncbi/SKESA
7. seqkit:https://github.com/shenwei356/seqkit/releases
-----
## Flags and Options
The pipeline provides several flgas to customize execution:
```sh
Usage: SGP.sh -g GENE_FILE -s SEROTYPE_FILE [--mode light|heavy] [--sra on|off] [--sra_number N] [-c COVERAGE] [-i IDENTITY] [-o OUTPUT_FILE]
-g GENE_FILE : FASTA file with target gene sequence (required).
-s SEROTYPE_FILE : File containing Salmonella serotypes (required).
--mode MODE : Choose between heavyweight mode and lightweight mode (default: light).
--sra on|off : Enable or disable SRA assembly (default: off).
--sra_num N : Specify the number of SRA files to retrieve (max: 100, default 50).
-c COVERAGE : The minimum genome coverage (default: 80%).
-i IDENTITY : The minimum percentage of identity (default: 90%).
-o OUTPUT_FILE : Output result file (default: gene_summary.tsv).
```
Example usage:
```sh
bash SGP.sh -g target_genes.fasta -s serotype_list.txt --mode heavy --sra on --sra_number 80 -c 95 -i 95 -o output.csv
```

## Output
Examples of the final output file:

Light mode
| Serotype            | Gene_ID                     | Min_coverage | Min_percentage_of_identity | Total_draft_genomes |  Total_complete_genomes | Complete_genomes_with_target_gene | Percentage_with_target_gene_in_complete_genomes |
|---------------------|-----------------------------|--------------|----------------------------|---------------------|-------------------------|-----------------------------------|-------------------------------------------------|
| Enteritidis         | NC_003197.2:1707344-1707789 | 95           | 95                         | 41408               | 165                     | 164                               | 99.00%                                          |
| Typhimurium         | NC_003197.2:1707344-1707789 | 95           | 95                         | 37382               | 226                     | 226                               | 100.00%                                         |
| Typhimurium var. 5- | NC_003197.2:1707344-1707789 | 95           | 95                         | 825                 | 11                      | 11                                | 100.00%                                         |
| Heidelberg          | NC_003197.2:1707344-1707789 | 95           | 95                         | 4905                | 51                      | 51                                | 100.00%                                         |
| Infantis            | NC_003197.2:1707344-1707789 | 95           | 95                         | 16915               | 80                      | 80                                | 100.00%                                         |
| Newport             | NC_003197.2:1707344-1707789 | 95           | 95                         | 8116                | 92                      | 92                                | 100.00%                                         |
| Uganda              | NC_003197.2:1707344-1707789 | 95           | 95                         | 1254                | 12                      | 12                                | 100.00%                                         |
| Braenderup          | NC_003197.2:1707344-1707789 | 95           | 95                         | 2565                | 5                       | 5                                 | 100.00%                                         |
| Muenchen            | NC_003197.2:1707344-1707789 | 95           | 95                         | 2384                | 22                      | 22                                | 100.00%                                         |
| Montevideo          | NC_003197.2:1707344-1707789 | 95           | 95                         | 4831                | 32                      | 31                                | 96.00%                                          |
| Javiana             | NC_003197.2:1707344-1707789 | 95           | 95                         | 1748                | 8                       | 8                                 | 100.00%                                         |
| Reading             | NC_003197.2:1707344-1707789 | 95           | 95                         | 1973                | 12                      | 12                                | 100.00%                                         |
| Dublin              | NC_003197.2:1707344-1707789 | 95           | 95                         | 3831                | 28                      | 28                                | 100.00%                                         |
| Oranienburg         | NC_003197.2:1707344-1707789 | 95           | 95                         | 1543                | 5                       | 5                                 | 100.00%                                         |
| Potsdam             | NC_003197.2:1707344-1707789 | 95           | 95                         | 128                 | 1                       | 1                                 | 100.00%                                         |
| Thompson            | NC_003197.2:1707344-1707789 | 95           | 95                         | 1727                | 17                      | 17                                | 100.00%                                         |
| Saintpaul           | NC_003197.2:1707344-1707789 | 95           | 95                         | 3873                | 24                      | 24                                | 100.00%                                         |
| Hadar               | NC_003197.2:1707344-1707789 | 95           | 95                         | 1927                | 15                      | 15                                | 100.00%                                         |
| Schwarzengrund      | NC_003197.2:1707344-1707789 | 95           | 95                         | 3038                | 27                      | 27                                | 100.00%                                         |
| Anatum              | NC_003197.2:1707344-1707789 | 95           | 95                         | 4843                | 40                      | 40                                | 100.00%                                         |
| Berta               | NC_003197.2:1707344-1707789 | 95           | 95                         | 604                 | 3                       | 3                                 | 100.00%                                        
