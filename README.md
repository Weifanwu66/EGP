# Gene Prevalence Estimation Tool for Enterobacteriaceae

## Overview
This tool is designed for estimating the prevalence of a specific gene in Enterobacteriaceae taxa, integrating NCBI genome retrieval, BLAST database construction, and automated query analysis.

## Features
**Genomic Data Acquisition**
### Genomic Data Acquisition
* **Storage‑friendly metadata:** Because of the large storage requirements, genome sequence files used to build the BLAST database are **not** stored in this repository. Instead, `database/metadata` holds the assembly‑accession lists so users can re‑download any sequence on demand.  
* **Pre‑built complete‑genome database:** A BLAST database built from complete genomes of the default seven Enterobacteriaceae genera is provided for quick, high‑quality searches.  
* **Heavy mode (`‑H heavy` + `‑t <taxon_file>`):** Adds draft genomes to the analysis. Draft assemblies for each target taxon are downloaded with *ncbi‑genome‑download*, shuffled, and sampled in iterations (Cochran’s formula with finite‑population correction; ≤ 20 iterations) to capture diversity while keeping runtime reasonable.  
* **Custom genome panel (`‑d <genus_file>`):** Lets users work **outside** the default Enterobacteriaceae set. Supply a text file (one genus per line) and the pipeline will download the corresponding **complete genomes**, build a bespoke BLAST database in real time, and then run either LIGHT or HEAVY mode against that database.

> **Note:** Building either the pre‑built archive or a large custom panel requires significant disk and CPU resources. Run these steps on an HPC system or a workstation with ≥ 250 GB free space.

---

**Genome files Organization**
Creates a structured directory hierarchy per genus → species → (for *Salmonella enterica*) subspecies → serotype. Unclassified genomes are placed in an `unclassified/` subfolder for each taxon.  
Complete‑genome downloads live under `complete_genome/`; draft‑genome iterations are stored separately under `draft_genome/`.  
The script to download complete genomes and their corresponding BLAST DB is build_EB_complete_genomes_database.sh
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
