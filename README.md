# Gene Prevalence Estimation Tool for Enterobacteriaceae

## Overview
This tool is designed for estimating the prevalence of a specific gene in Enterobacteriaceae taxa, integrating NCBI genome retrieval, BLAST database construction, and automated query analysis.

## Features
**Genomic Data Acquisition**
  - Due to the large storage requirements, genome sequence files used to build BLAST database will not be uploaded to this repository. Instead, metadata files containing all assembly accessions for the downloaded genomes are provided for each respective directory and stored in `database/metadata`. This ensures traceability and allows users to retrieve specific assemblies if needed.
  - A pre-built BLAST database has been constructed for complete genomes. However, since some complete genomes are relatively small and may not be representative of the full genetic diversity of a taxon, users may choose to enable **heavy mode** to include draft genomes in their analysis.
  - For **heavy mode**, users must provide their target taxa. Draft genomes for each taxon are retrieved and randomly shuffled before selection based on their accessions using ncbi-genome-download. The draft genomes are iteratively sampled due to their large number, ensuring representative sampling across taxa.
  - The sample size per iteration is automatically calculated using Cochran’s formula and finite population correction and and number of iterations are determined using square-root scaling, based on the total number of draft genomes (contigs) available in GenBank. To maintain computational feasibility, the number of iterations is capped at 20.

**Genome files Organization**
- Creates structured directories per genus, species, *Salmonella enterica* subspecies, and serotypes under *Salmonella enterica subsp. enterica*.
- Maintains an unclassified directory for each taxon.
- Draft genomes are downloaded randomly in iterations and stored separately within the draft genome directory.
- The script to download complete genomes is download_complete_genomes.sh and the script to build the BLAST database is makeblastdb_complete_genomes.sh.
```
│   ├── Escherichia/
│   │   ├── unclassified/
│   │   ├── Escherichia_coli/
│   │   ├── Escherichia_fergusonii/
│   │   ├── Escherichia_albertii/
│   │   ├── ...
│   ├── Salmonella/
│   │   ├── unclassified/
│   │   ├── Salmonella_enterica/
│   │   │   ├── unclassified/
│   │   │   ├── subsp_enterica/
│   │   │   │   ├── unclassified/
│   │   │   │   ├── Typhimurium/
│   │   │   │   ├── Infantis/
│   │   │   │   ├── Newport/
│   │   │   │   ├── Heidelberg/
│   │   │   │   ├── ...
│   │   │   ├── subsp_salamae/
│   │   │   ├── subsp_arizonae/
│   │   │   ├── subsp_diarizonae/
│   │   │   ├── subsp_houtenae/
│   │   │   ├── subsp_indica/
│   │   ├── Salmonella_bongori/
│   ├── Shigella/
│   │   ├── unclassified/
│   │   ├── Shigella_flexneri/
│   │   ├── Shigella_sonnei/
│   │   ├── Shigella_boydii/
│   │   ├── ...
│   ├── Klebsiella/
│   │   ├── unclassified/
│   │   ├── Klebsiella_pneumoniae/
│   │   ├── Klebsiella_oxytoca/
│   │   ├── ...
│   ├── Enterobacter/
│   │   ├── unclassified/
│   │   ├── Enterobacter_cloacae/
│   │   ├── Enterobacter_hormaechei/
│   │   ├── ...
│   ├── Citrobacter/
│   │   ├── unclassified/
│   │   ├── Citrobacter_freundii/
│   │   ├── Citrobacter_koseri/
│   │   ├── ...
│   ├── Cronobacter/
│   │   ├── unclassified/
│   │   ├── Cronobacter_sakazakii/
│   │   ├── Cronobacter_malonaticus/
│   │   ├── ...
│   ├── Proteus/
│   │   ├── complete_genomes/
│   │   ├── Proteus_mirabilis/
│   │   ├── Proteus_vulgaris/
│   │   ├── ...
```
**BLAST Query & Analysis**
- BLAST analysis can be run in two modes:
  1. **Light mode**: The analysis is limited to pre-built complete genome databases, leveraging high-quality genome assemblies for faster performance.
  2. **Heavy mode**: The pipeline runs BLAST searches against both complete and draft genome databases. The draft genomes databases are dynamically constructed during runtime to ensure representative sampling across each target group.
- Query gene file is provided by users (supports batch processing of multiple genes).
- Filters results by user-defined minimum identity & coverage thresholds.

**Gene prevalence calculation**
- **Light mode**: The prevalence of the gene is calculated based on hits in **complete genomes only**.
- **heavy mode**: In addition to analyzing **complete genomes**, this mode incorporates **Draft genomes** to provide a more comprehensive prevalence estimate. Since draft genomes may have variable quality, multiple iterations may perform. The final prevalence estimate is averaged over all iterations, ensuring robustness against sampling bias.
------
## Installation
To run this pipeline, set up a Conda environment with the required dependencies.
1. Clone the Repository
```sh
git clone https://github.com/Weifanwu66/EGP.git
cd EGP
```
2. Create and Activate the Conda Environment
The pipeline requires a Conda environment with all necessary dependencies. To create and activate it, run:
```sh
conda env create -f environment.yml
conda activate EGP
```
To verify the installation, check if all tools are installed:
```sh
conda list | grep -E "blast|ncbi-genome-download|entrez-direct"
```
If any package is missing, please install it manually:
```sh
conda install -c bioconda <package_name>
```
-----
## Dependencies
1. ncbi-genome-download: Blin, K. (2023). ncbi-genome-download (0.3.3). Zenodo. https://doi.org/10.5281/zenodo.8192486
2. NCBI BLAST: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html
3. entrez-direct: https://www.ncbi.nlm.nih.gov/books/NBK179288/
-----
## Example usage:
To download in default light mode
```sh
bash EGP.sh -g target_genes.fasta -t taxa_list.txt
```
If no taxa is provided, it will automatically screen for the whole database
```sh
bash EGP.sh -g target_genes.fasta
```
To define your minimum coverage and identity
```sh
bash EGP.sh -g target_genes.fasta -t taxon_list.txt -c 95 -i 95
```
To turn on heavy mode
```sh
bash EGP.sh -g target_genes.fasta -t taxon_list.txt --mode heavy
```
If you just want to look at one organism
```sh
bash EGP.sh -g target_genes.fasta -t "Salmonella"
```
You can choose to turn on heavy mode and name your output file, but by default, the output file will be stored in `result/gene_summary.csv`.
```sh
bash EGP.sh -g target_genes.fasta -t taxon_list.txt --mode heavy -c 95 -i 95 -o output.csv
```
To overwrite the previous result, add `--overwrite`
## Output
Examples of the final output file:

**Light mode**
| Organism                            | Gene_ID                      | Min_percentage_of_coverage | Min_percentage_of_identity | Total_draft_genomes | Total_complete_genomes | Complete_genomes_with_target_genes | Percentage_with_target_genes_complete_genomes |
| ----------------------------------- | ---------------------------- | -------------------------- | -------------------------- | ------------------- | ---------------------- | ---------------------------------- | --------------------------------------------- |
| Salmonella enterica subsp. enterica | NC_003197.2:c3010816-3009905 | 90                         | 95                         | 235792              | 2307                   | 1924                               | 83.00%                                        |
| Salmonella enterica subsp. enterica | NC_003197.2:3019846-3021524  | 90                         | 95                         | 235792              | 2307                   | 2270                               | 98.00%                                        |
| Salmonella enterica subsp. enterica | NC_003197.2:c2925778-2923593 | 90                         | 95                         | 235792              | 2307                   | 2276                               | 98.00%                                        |

**Heavy mode**
| Organism          | Gene_ID                     | Min_percentage_of_coverage | Min_percentage_of_identity | Total_draft_genomes | Total_complete_genomes | Draft_genomes_sample_size | Number_of_iterations | Complete_genomes_with_target_genes | Draft_genomes_with_target_genes | Percentage_with_target_genes_complete_genomes | Percentage_with_target_genes_draft_genomes |
| ----------------- | --------------------------- | -------------------------- | -------------------------- | ------------------- | ---------------------- | ------------------------- | -------------------- | ---------------------------------- | ------------------------------- | --------------------------------------------- | ------------------------------------------ |
| Salmonella Uganda | NC_003197.2:1707344-1707789 | 80                         | 90                         | 1173                | 12                     | 290                       | 1                    | 12                                 | 290                             | 100.00%                                       | 100.00%                                    |
