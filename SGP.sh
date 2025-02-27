#!/bin/bash
# set default values
MODE="light"
MIN_COVERAGE=80
MIN_IDENTITY=90
SRA_FLAG="off"
SRA_NUM=50
SEROTYPE_FILE=""
GENE_FILE=""
OUTPUT_FILE="gene_summary.csv"
ALL_GENOMES_FILE="all_genomes.tsv"
ASSEMBLY_DIR="assemblies"
# usage setup
usage(){
echo "Usage: $0 -g GENE_FILE -s SEROTYPE_FILE [--mode light|heavy] [--sra on|off] [--sra_number N] [-c COVERAGE] [-i IDENTITY] [-o OUTPUT_FILE]"
echo "-g GENE_FILE : FASTA file with target gene sequence (required)." 
echo "-s SEROTYPE_FILE : File containing Salmonella serotypes (required)."
echo "--mode MODE : Choose between heavyweight mode and lightweight mode (default: lightweight)."
echo "--sra on|off : Enable or disable SRA assembly (default: off)."
echo "--sra_num N : Specify the number of SRA files to retrieve (max: 100, default 50)."
echo "-c COVERAGE : The minimum genome coverage (default: 80%)." 
echo "-i IDENTITY : The minimum percentage of identity (default: 90%)."
echo "-o OUTPUT_FILE : Output result file (default: gene_summary.tsv)."
}

# Parse argument
while [[ "$#" -gt 0 ]]; do
case "$1" in
-s) SEROTYPE_FILE="$2"; shift ;;
-g) GENE_FILE="$2"; shift ;;
-c) MIN_COVERAGE="$2"; shift ;;
-i) MIN_IDENTITY="$2"; shift ;;
-o) OUTPUT_FILE="$2"; shift ;;
--mode) MODE="$2"; shift ;;
--sra) SRA_FLAG="$2"; shift ;;
--sra_number) SRA_NUM="$2"
if [[ "$SRA_NUM" -gt 100 ]]; then
echo "Maximum allowed SRA file is 100. Setting --sra_number to 100."
SRA_NUM=100
fi
shift ;;
*) echo "Invalid option" usage; exit 1 ;;
esac
shift
done

#Ensure both serotype list file and gene files are provided
if [[ -z "$GENE_FILE" || -z "$SEROTYPE_FILE" ]]; then
echo "Error! Please provide a gene sequence file and your serotype list!"
usage
exit 1
fi
# Output directory setup
mkdir -p genomes blast_results "$ASSEMBLY_DIR"
# Retrieve assembly accessions based on serotypes
echo -e "Serotype\tAccessions" > "$ALL_GENOMES_FILE"
mapfile -t SEROTYPES < "$SEROTYPE_FILE"
for SEROTYPE in "${SEROTYPES[@]}"; do
CLEAN_SEROTYPE=$(echo "$SEROTYPE" | sed 's/\.//g; s/,//g; s/ /_/g; s/-/_/g; s/\[//g; s/\]//g')
echo "Retrieving assembly accessions for $SEROTYPE"
mkdir -p "genomes/${CLEAN_SEROTYPE}"
esearch -db assembly -query "Salmonella enterica AND $SEROTYPE" | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession,AssemblyStatus,Organism > "genomes/${CLEAN_SEROTYPE}_accessions_raw.txt"
# Filter out Typhimurium var. 5- from Typhimurium
if [[ "$SEROTYPE" == "Typhimurium" ]]; then
grep -v "var. 5-" "genomes/${CLEAN_SEROTYPE}_accessions_raw.txt" > "genomes/${CLEAN_SEROTYPE}_accessions.txt"
else
mv "genomes/${CLEAN_SEROTYPE}_accessions_raw.txt" "genomes/${CLEAN_SEROTYPE}_accessions.txt"
fi
# Define accessions for complete genomes and save them to all_genomes.tsv"
grep "Complete Genome" "genomes/${CLEAN_SEROTYPE}_accessions.txt" | awk -v serotype="$CLEAN_SEROTYPE" '{print serotype "\t" $1}' >> "$ALL_GENOMES_FILE"
grep -v "Complete Genome" "genomes/${CLEAN_SEROTYPE}_accessions.txt" > "genomes/${CLEAN_SEROTYPE}_draft_genomes.txt"
rm "genomes/${CLEAN_SEROTYPE}_accessions.txt"
done
# Download genome files based on accessions per serotype
mapfile -t ACCESSIONS < <(awk 'NR>1' "$ALL_GENOMES_FILE")
for LINE in "${ACCESSIONS[@]}"; do
IFS=$'\t' read -r CLEAN_SEROTYPE ACCESSION <<< "$LINE"
ZIP_FILE="genomes/$CLEAN_SEROTYPE/$ACCESSION.zip"
UNZIP_DIR="genomes/$CLEAN_SEROTYPE"
if ls "$UNZIP_DIR/${ACCESSION}"*"_genomic.fna" 1> /dev/null 2>&1; then
echo "Genome already exists for $CLEAN_SEROTYPE ($ACCESSION), skipping downloading"
continue
fi
echo "Downloading genomes for $CLEAN_SEROTYPE"
datasets download genome accession "$ACCESSION" --no-progressbar --filename "$ZIP_FILE" &> /dev/null < /dev/null
echo "Unzipping genome for $CLEAN_SEROTYPE"
unzip -q -o "$ZIP_FILE" -d "$UNZIP_DIR"
mv "$UNZIP_DIR/ncbi_dataset/data/"*/*_genomic.fna "$UNZIP_DIR/" 2>/dev/null
rm "$ZIP_FILE"
rm -rf "$UNZIP_DIR/ncbi_dataset"
done 
# Set up the output file headers
if [[ "$MODE" == "heavy" ]]; then
echo -e "Serotype,Gene_ID,Min_percentage_of_coverage,Min_percentage_of_identity,Total_draft_genomes,Total_complete_genomes,SRA_requested,SRA_dropped_low_assemble_quality,Total_assembled_genomes,Complete_genomes_with_target_genes,Assembled_genomes_with_target_genes,Percentage_with_target_genes_complete_genomes,Percentage_with_target_genes_assembled_genomes" > "$OUTPUT_FILE"
else
echo -e "Serotype,Gene_ID,Min_percentage_of_coverage,Min_percentage_of_identity,Total_draft_genomes,Total_complete_genomes,Complete_genomes_with_target_genes,Percentage_with_target_genes_complete_genomes" > "$OUTPUT_FILE"
fi
# Run BLAST on all genomes to check gene presence per serotype
mapfile -t SEROTYPES < "$SEROTYPE_FILE"
for SEROTYPE in "${SEROTYPES[@]}"; do
CLEAN_SEROTYPE=$(echo "$SEROTYPE" | sed 's/ /_/g;s/\.//g;s/-/_/g;s/,//g')
GENOME_DIR="genomes/$CLEAN_SEROTYPE"
mkdir -p "$GENOME_DIR"
CONCATENATED_GENOMES="genomes/all_${CLEAN_SEROTYPE}.fna"
if [[ -f "$CONCATENATED_GENOMES" && -s "$CONCATENATED_GENOMES" && -f "blast_results/${CLEAN_SEROTYPE}_db.nsq" ]]; then
echo "BLAST database already exists for $CLEAN_SEROTYPE, skipping"
else
echo "Making BLAST database and running BLAST for $SEROTYPE"
find "$GENOME_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$CONCATENATED_GENOMES"
makeblastdb -in "$CONCATENATED_GENOMES" -dbtype nucl -out "blast_results/${CLEAN_SEROTYPE}_db"
blastn -query "$GENE_FILE" -db "blast_results/${CLEAN_SEROTYPE}_db" -out "blast_results/${CLEAN_SEROTYPE}_results.tsv" -outfmt 6 -perc_identity "$MIN_IDENTITY"
fi
awk -v min_cov="$MIN_COVERAGE" '(($4/($8 - $7 + 1)) * 100) >= min_cov' "blast_results/${CLEAN_SEROTYPE}_results.tsv" > "blast_results/filtered_${CLEAN_SEROTYPE}_results.tsv"
TOTAL_DRAFT_GENOMES=$(awk 'NR>1' "genomes/${CLEAN_SEROTYPE}_draft_genomes.txt" 2>/dev/null | wc -l)
TOTAL_COMPLETE_GENOMES=$(find "$GENOME_DIR" -type f -name "*_genomic.fna" | wc -l)
# When heavy mode is on and complete genomes are less than 30
if [[ "$MODE" == "heavy" && "$TOTAL_COMPLETE_GENOMES" -lt 30 && "$SRA_FLAG" == "on" ]]; then
SRA_DIR="$ASSEMBLY_DIR/${CLEAN_SEROTYPE}"
mkdir -p "$SRA_DIR"
if [[ -f "$SRA_DIR/sra_accessions.txt" ]]; then
echo "Skipping SRA accession retrieval; using existing $SRA_DIR/sra_accessions.txt"
else
echo "Processing assemblies for serotypes with total complete genomes under 30, now processing $SEROTYPE"
esearch -db sra -query "Salmonella enterica AND ${SEROTYPE}" | efetch -format runinfo | awk -F"," 'NR>1 {print $1}' | head -n "$SRA_NUM" > "$SRA_DIR/sra_accessions.txt"
fi
mapfile -t SRA_ACCESSIONS < "$SRA_DIR/sra_accessions.txt"
for SRA_ACC in "${SRA_ACCESSIONS[@]}"; do
# Skip processing if all expected FASTQ files exist
if [[ -f "$SRA_DIR/${SRA_ACC}_paired_1.fastq" && -f "$SRA_DIR/${SRA_ACC}_paired_2.fastq" && \
      -f "$SRA_DIR/${SRA_ACC}_unpaired_1.fastq" && -f "$SRA_DIR/${SRA_ACC}_unpaired_2.fastq" ]]; then
     echo "Skipping assembly for $SRA_ACC; FASTQ files already exist."
continue
fi
prefetch "$SRA_ACC" --output-directory "$SRA_DIR"
fasterq-dump "$SRA_DIR/$SRA_ACC/$SRA_ACC.sra" --outdir "$SRA_DIR" --split-files --threads 16
SEQ_PLATFORM=$(esearch -db sra -query "$SRA_ACC" | efetch -format runinfo | awk -F"," 'NR==2 {print $20}')
# Map sequencing platform to adapters
declare -A SRA_ADAPTER_MAP
SRA_ADAPTER_MAP["Illumina MiSeq"]="TruSeq3-PE.fa"
SRA_ADAPTER_MAP["Illumina HiSeq"]="TruSeq3-PE.fa"
SRA_ADAPTER_MAP["Illumina NextSeq"]="NexteraPE-PE.fa"
DEFAULT_ADAPTER="TruSeq3-PE.fa"
ADAPTERS="${SRA_ADAPTER_MAP[$SEQ_PLATFORM]:-$DEFAULT_ADAPTER}"
echo "Running trimmomatic"
 trimmomatic PE  \
"$SRA_DIR/${SRA_ACC}_1.fastq" "$SRA_DIR/${SRA_ACC}_2.fastq" \
"$SRA_DIR/${SRA_ACC}_paired_1.fastq" "$SRA_DIR/${SRA_ACC}_unpaired_1.fastq" \
"$SRA_DIR/${SRA_ACC}_paired_2.fastq" "$SRA_DIR/${SRA_ACC}_unpaired_2.fastq" \
ILLUMINACLIP:"$ADAPTERS":2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
echo "Assembling with Skesa"
skesa --reads "$SRA_DIR/${SRA_ACC}_paired_1.fastq,$SRA_DIR/${SRA_ACC}_paired_2.fastq" --contigs_out "$SRA_DIR/$SRA_ACC/contigs.fasta" --use_paired_ends --min_contig 500 
# Quality check by using seqkit
echo "Evaluating assembly quality for $SRA_ACC"
SEQKIT_STATS=$(seqkit stats --all --tabular "$SRA_DIR/$SRA_ACC/contigs.fasta" | awk 'NR==2')
N50=$(echo "$SEQKIT_STATS" | awk '{print $13}')
TOTAL_CONTIGS=$(echo "$SEQKIT_STATS" | awk '{print $4}')
LARGEST_CONTIG=$(echo "$SEQKIT_STATS" | awk '{print $8}')
if [[ "$N50" -ge 50000 && "$TOTAL_CONTIGS" -lt 100 && "$LARGEST_CONTIG" -ge 200000 ]]; then
echo "Assembly $SRA_ACC passed quality control"
else
echo "Assembly $SRA_ACC failed quality control, discarding"
rm -rf "$SRA_DIR/$SRA_ACC" #removed the failed assembly
continue
fi
echo "Running blast for assembled genome for $SEROTYPE"
awk -v sra_acc="$SRA_ACC" '/^>/{$0=">" sra_acc "_" substr($0, 2)}1' "$SRA_DIR/$SRA_ACC/contigs.fasta" > "$SRA_DIR/$SRA_ACC/modified_contigs.fasta"
mv "$SRA_DIR/$SRA_ACC/modified_contigs.fasta" "$SRA_DIR/$SRA_ACC/contigs.fasta"
cat "$SRA_DIR/$SRA_ACC/contigs.fasta" >> "$SRA_DIR/combined_assembled_genomes"
makeblastdb -in "$SRA_DIR/combined_assembled_genomes" -dbtype nucl -out "blast_results/${CLEAN_SEROTYPE}_assembled_db"
blastn -query "$GENE_FILE" -db "blast_results/${CLEAN_SEROTYPE}_assembled_db" -out "blast_results/${CLEAN_SEROTYPE}_assembled_results.tsv" -outfmt 6 -perc_identity "$MIN_IDENTITY"
awk -v min_cov="$MIN_COVERAGE" '(($4/($8 - $7 + 1))*100) >= min_cov' "blast_results/${CLEAN_SEROTYPE}_assembled_results.tsv" > "blast_results/filtered_${CLEAN_SEROTYPE}_assembled_results.tsv"
done
fi
done
# Ensure serotypes are loaded and process results for each serotype
mapfile -t SEROTYPES < "$SEROTYPE_FILE"
for SEROTYPE in "${SEROTYPES[@]}"; do
CLEAN_SEROTYPE=$(echo "$SEROTYPE" | sed 's/ /_/g; s/-/_/g; s/\.//g; s/,//g')
TOTAL_COMPLETE_GENOMES=$(find "genomes/$CLEAN_SEROTYPE" -type f -name "*_genomic.fna" | wc -l)
TOTAL_DRAFT_GENOMES=$(awk 'NR>1' "genomes/${CLEAN_SEROTYPE}_draft_genomes.txt" 2>/dev/null | wc -l)
if [[ "$MODE" == "heavy" && "$TOTAL_COMPLETE_GENOMES" -lt 30 && "$SRA_FLAG" == "on" ]]; then
SRA_REQUESTED="$SRA_NUM"
TOTAL_ASSEMBLED_GENOMES=$(find "$ASSEMBLY_DIR/$CLEAN_SEROTYPE" -mindepth 1 -type d | wc -l)
SRA_DROPPED=$((SRA_REQUESTED - TOTAL_ASSEMBLED_GENOMES))
else
SRA_REQUESTED=0
TOTAL_ASSEMBLED_GENOMES=0
SRA_DROPPED=0
fi
# Extract all unique gene IDs detected either in complete or assembled genomes
ALL_GENES=$(grep "^>" "$GENE_FILE" | sed 's/>//' | awk '{print $1}')
GENE_WITH_HITS=$(awk '{print $1}' "blast_results/filtered_${CLEAN_SEROTYPE}_results.tsv" 2>/dev/null; awk '{print $1}' "blast_results/filtered_${CLEAN_SEROTYPE}_assembled_results.tsv" 2>/dev/null | sort -u)
for GENE_ID in $ALL_GENES; do
if grep -q "^${GENE_ID}$" <<< "$GENE_WITH_HITS"; then
COMPLETE_GENOMES_WITH_TARGET_GENES=$(awk -v gene="$GENE_ID" '$1 == gene {print $2}' "blast_results/filtered_${CLEAN_SEROTYPE}_results.tsv" 2>/dev/null | sort -u | wc -l)
if [[ "$MODE" == "heavy" && "$TOTAL_COMPLETE_GENOMES" -lt 30 && "$SRA_FLAG" == "on" ]]; then
ASSEMBLED_GENOMES_WITH_TARGET_GENES=$(awk -v gene="$GENE_ID" '$1 == gene {print $2}' "blast_results/filtered_${CLEAN_SEROTYPE}_assembled_results.tsv" 2>/dev/null | awk -F'_' '{print $1}' | sort -u | wc -l)
else
ASSEMBLED_GENOMES_WITH_TARGET_GENES=0
fi
if [[ "$TOTAL_ASSEMBLED_GENOMES" -gt 0 ]]; then
PERCENT_WITH_TARGET_GENES_ASSEMBLED_GENOMES=$(echo "scale=2; ($ASSEMBLED_GENOMES_WITH_TARGET_GENES/$TOTAL_ASSEMBLED_GENOMES) * 100" | bc 2>/dev/null)
else
PERCENT_WITH_TARGET_GENES_ASSEMBLED_GENOMES=0
fi
if [[ "$TOTAL_COMPLETE_GENOMES" -gt 0 ]]; then
PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES=$(echo "scale=2; ($COMPLETE_GENOMES_WITH_TARGET_GENES/$TOTAL_COMPLETE_GENOMES) * 100" | bc 2>/dev/null)
else
PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES=0
fi
if [[ "$MODE" == "heavy" ]]; then
echo -e "$SEROTYPE,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$SRA_REQUESTED,$SRA_DROPPED,$TOTAL_ASSEMBLED_GENOMES,$COMPLETE_GENOMES_WITH_TARGET_GENES,$ASSEMBLED_GENOMES_WITH_TARGET_GENES,$PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES%,$PERCENT_WITH_TARGET_GENES_ASSEMBLED_GENOMES%" >> "$OUTPUT_FILE"
else
echo -e "$SEROTYPE,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$COMPLETE_GENOMES_WITH_TARGET_GENES,$PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES%"  >> "$OUTPUT_FILE"
fi
else
# If no hits were found for the gene, output 0 %
if [[ "$MODE" == "heavy" ]]; then
echo -e "$SEROTYPE,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$SRA_REQUESTED,$SRA_DROPPED,$TOTAL_ASSEMBLED_GENOMES,0,0,0%,0%" >> "$OUTPUT_FILE"
else
echo -e "$SEROTYPE,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,0,0%" >> "$OUTPUT_FILE"
fi
fi
done
done
echo "Analysis complete. Results saved in $OUTPUT_FILE"
