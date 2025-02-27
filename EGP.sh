#!/bin/bash
# set default values
MODE="light"
MIN_COVERAGE=80
MIN_IDENTITY=90
SRA_FLAG="off"
SRA_NUM=50
SPECIES_SEROTYPE_FILE=""
GENE_FILE=""
OUTPUT_FILE="EB_gene_summary.csv"
ALL_GENOMES_FILE="EB_all_genomes.tsv"
ASSEMBLY_DIR="EB_assemblies"
# usage setup
usage(){
echo "Usage: $0 -g GENE_FILE -s SEROTYPE_FILE [--mode light|heavy] [--sra on|off] [--sra_number N] [-c COVERAGE] [-i IDENTITY] [-o OUTPUT_FILE]"
echo "-g GENE_FILE : FASTA file with target gene sequence (required)." 
echo "-s SPECIES_SEROTYPE_FILE : File containing species and serotypes if species is Salmonella enterica (required)."
echo "--mode MODE : Choose between heavyweight mode and lightweight mode (default: light)."
echo "--sra on|off : Enable or disable SRA assembly (default: off)."
echo "--sra_num N : Specify the number of SRA files to retrieve (max: 100, default 50)."
echo "-c COVERAGE : The minimum genome coverage (default: 80%)." 
echo "-i IDENTITY : The minimum percentage of identity (default: 90%)."
echo "-o OUTPUT_FILE : Output result file (default: EB_gene_summary.tsv)."
}

# Parse argument
while [[ "$#" -gt 0 ]]; do
case "$1" in
-s) SPECIES_SEROTYPE_FILE="$2"; shift ;;
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
if [[ -z "$GENE_FILE" || -z "$SPECIES_SEROTYPE_FILE" ]]; then
echo "Error! Please provide a gene sequence file and your species/serotype list!"
usage
exit 1
fi
# Output directory setup
mkdir -p EB_genomes EB_blast_results "$ASSEMBLY_DIR"
# Retrieve assembly accessions based on serotypes
echo -e "Species\tSerotype\tAccessions" > "$ALL_GENOMES_FILE"
sed -i 's/ \+/\t/g' "$SPECIES_SEROTYPE_FILE"
awk '{
species = $1 " " $2;
serotype = (NF > 2) ? substr($0, index($0, $3)) : "NA";
print species "\t" serotype;
}' "$SPECIES_SEROTYPE_FILE" > temp_file && mv temp_file "$SPECIES_SEROTYPE_FILE"
mapfile -t LINES < "$SPECIES_SEROTYPE_FILE"
for LINE in "${LINES[@]}"; do
IFS=$'\t' read -r SPECIES SEROTYPE <<< "$LINE"
SPECIES=$(echo "$SPECIES" | xargs)
SEROTYPE=$(echo "$SEROTYPE" | xargs)
CLEAN_SEROTYPE=$(echo "$SEROTYPE" | sed 's/\.//g; s/,//g; s/ /_/g; s/-/_/g; s/\[//g; s/\]//g')
CLEAN_SPECIES=$(echo "$SPECIES" | sed 's/\.//g; s/,//g; s/ /_/g; s/-/_/g; s/\[//g; s/\]//g')
if [[ "$SEROTYPE" == "NA" ]]; then
CLEAN_SEROTYPE=""
echo "Downloading complete genomes for $SPECIES directly using ncbi-datasets tool"
GENOME_DIR="EB_genomes/${CLEAN_SPECIES}"
mkdir -p "$GENOME_DIR"
datasets download genome taxon "$SPECIES" --assembly-level complete --filename "EB_genomes/${CLEAN_SPECIES}.zip"
if [[ ! -s "EB_genomes/${CLEAN_SPECIES}.zip" ]]; then
echo "ERROR: No genomes downloaded for $SPECIES"
continue
fi
unzip -q -o "EB_genomes/${CLEAN_SPECIES}.zip" -d "$GENOME_DIR"
mv "${GENOME_DIR}/ncbi_dataset/data/"*/*_genomic.fna "${GENOME_DIR}" 2>/dev/null
rm -rf "${GENOME_DIR}/ncbi_dataset"
rm -rf "EB_genomes/${CLEAN_SPECIES}.zip"
else
echo "Retrieving assembly accessions for $SPECIES $SEROTYPE using esearch"
QUERY="$SPECIES[Organism] AND $SEROTYPE[All Fields]"
GENOME_DIR="EB_genomes/${CLEAN_SPECIES}_${CLEAN_SEROTYPE}"
mkdir -p "$GENOME_DIR"
esearch -db assembly -query "$QUERY AND complete genome [filter]" | efetch -format docsum | xtract -pattern DocumentSummary -element AssemblyAccession,Organism > "${GENOME_DIR}/accessions.txt"
# Filter out Typhimurium var. 5- from Typhimurium
if [[ "$SPECIES" == "Salmonella enterica" && "$SEROTYPE" == "Typhimurium" ]]; then
grep -v "var. 5-" "${GENOME_DIR}/accessions.txt" > "${GENOME_DIR}/filtered_accessions.txt"
mv "${GENOME_DIR}/filtered_accessions.txt" "${GENOME_DIR}/accessions.txt"
fi
# Store accessions with their corresponding organisms in all_genomes.tsv file
cat "${GENOME_DIR}/accessions.txt" | awk -v species="$CLEAN_SPECIES" -v serotype="$CLEAN_SEROTYPE" '{print species"\t"serotype"\t"$1}' >> "$ALL_GENOMES_FILE"
# Download genome files based on accessions
mapfile -t LINES < "$ALL_GENOMES_FILE"
for LINE in "${LINES[@]}"; do
IFS=$'\t' read -r CLEAN_SPECIES CLEAN_SEROTYPE ACCESSION <<< "$LINE"
GENOME_DIR="EB_genomes/${CLEAN_SPECIES}_${CLEAN_SEROTYPE}"
ZIP_FILE="$GENOME_DIR/$ACCESSION.zip"
if ls "$GENOME_DIR/${ACCESSION}"*"_genomic.fna" 1> /dev/null 2>&1; then
echo "Genome already exists for $CLEAN_SPECIES $CLEAN_SEROTYPE ($ACCESSION), skipping downloading"
continue
fi
echo "Downloading genomes for $CLEAN_SPECIES $CLEAN_SEROTYPE"
datasets download genome accession "$ACCESSION" --no-progressbar --filename "$ZIP_FILE" &> /dev/null < /dev/null
unzip -q -o "$ZIP_FILE" -d "$GENOME_DIR"
mv "$GENOME_DIR/ncbi_dataset/data/"*/*_genomic.fna "$GENOME_DIR" 2>/dev/null
rm "$ZIP_FILE" 
rm -rf "$GENOME_DIR/ncbi_dataset"
done
fi
done
# Set up the output file headers
if [[ "$MODE" == "heavy" ]]; then
echo -e "Species,Serotype,Gene_ID,Min_percentage_of_coverage,Min_percentage_of_identity,Total_draft_genomes,Total_complete_genomes,SRA_requested,SRA_dropped,Total_assembled_genomes,Complete_genomes_with_target_genes,Assembled_genomes_with_target_genes,Percentage_with_target_genes_complete_genomes,Percentage_with_target_genes_assembled_genomes" > "$OUTPUT_FILE"
else
echo -e "Species,Serotype,Gene_ID,Min_percentage_of_coverage,Min_percentage_of_identity,Total_draft_genomes,Total_complete_genomes,Complete_genomes_with_target_genes,Percentage_with_target_genes_complete_genomes" > "$OUTPUT_FILE"
fi
# Run BLAST on all genomes to check gene presence per serotype
mapfile -t LINES < "$SPECIES_SEROTYPE_FILE"
for LINE in "${LINES[@]}"; do
IFS=$'\t' read -r SPECIES SEROTYPE <<< "$LINE"
CLEAN_SPECIES=$(echo "$SPECIES" | sed 's/ /_/g')
CLEAN_SEROTYPE=$(echo "$SEROTYPE" | sed 's/ /_/g;s/\.//g;s/-/_/g;s/,//g')
if [[ "$SEROTYPE" == "NA" ]]; then
CLEAN_SEROTYPE=""
GENOME_DIR="EB_genomes/${CLEAN_SPECIES}"
else
GENOME_DIR="EB_genomes/${CLEAN_SPECIES}_${CLEAN_SEROTYPE}"
fi
mkdir -p "$GENOME_DIR"
CONCATENATED_GENOMES="$GENOME_DIR/combined.fna"
find "$GENOME_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$CONCATENATED_GENOMES"
echo "Making blast database and running blast for $SPECIES $SEROTYPE"
if [[  -n "$CLEAN_SEROTYPE" ]]; then
BLAST_DB="EB_blast_results/${CLEAN_SPECIES}_${CLEAN_SEROTYPE}_db"
BLAST_OUTPUT="EB_blast_results/${CLEAN_SPECIES}_${CLEAN_SEROTYPE}_results.tsv"
FILTERED_BLAST_OUTPUT="EB_blast_results/filtered_${CLEAN_SPECIES}_${CLEAN_SEROTYPE}_results.tsv"
else
BLAST_DB="EB_blast_results/${CLEAN_SPECIES}_db"
BLAST_OUTPUT="EB_blast_results/${CLEAN_SPECIES}_results.tsv"
FILTERED_BLAST_OUTPUT="EB_blast_results/filtered_${CLEAN_SPECIES}_results.tsv"
fi
makeblastdb -in "$CONCATENATED_GENOMES" -dbtype nucl -out "$BLAST_DB"
blastn -query "$GENE_FILE" -db "$BLAST_DB" -out "$BLAST_OUTPUT" -outfmt 6 -perc_identity "$MIN_IDENTITY"
sync
awk -v min_cov="$MIN_COVERAGE" '(($4/($8 - $7 + 1)) * 100) >= min_cov' "$BLAST_OUTPUT" > "$FILTERED_BLAST_OUTPUT" 
TOTAL_COMPLETE_GENOMES=$(find "$GENOME_DIR" -type f -name "*_genomic.fna" | wc -l)
# When heavy mode is on and complete genomes are less than 30
if [[ "$MODE" == "heavy" && "$TOTAL_COMPLETE_GENOMES" -lt 30 && "$SRA_FLAG" == "on" ]]; then
SRA_DIR="$ASSEMBLY_DIR/${CLEAN_SPECIES}_${CLEAN_SEROTYPE}"
mkdir -p "$SRA_DIR"
echo "Processing assemblies for Salmonella serotypes with total complete genomes under 30, now processing $SEROTYPE"
esearch -db sra -query "${SPECIES} AND ${SEROTYPE}" | efetch -format runinfo | awk -F"," 'NR>1 {print $1}' | head -n "$SRA_NUM" > "$SRA_DIR/sra_accessions.txt"
mapfile -t SRA_ACCESSIONS < "$SRA_DIR/sra_accessions.txt"
for SRA_ACC in "${SRA_ACCESSIONS[@]}"; do
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
cat "$SRA_DIR/$SRA_ACC/contigs.fasta" >> "$SRA_DIR/combined_assembled_genomes"
makeblastdb -in "$SRA_DIR/combined_assembled_genomes" -dbtype nucl -out "EB_blast_results/${CLEAN_SPECIES}_${CLEAN_SEROTYPE}_assembled_db"
blastn -query "$GENE_FILE" -db "EB_blast_results/${CLEAN_SPECIES}_${CLEAN_SEROTYPE}_assembled_db" -out "EB_blast_results/${CLEAN_SPECIES}_${CLEAN_SEROTYPE}_assembled_results.tsv" -outfmt 6 -perc_identity "$MIN_IDENTITY"
awk -v min_cov="$MIN_COVERAGE" '(($4/($8 - $7 + 1))*100) >= min_cov' "EB_blast_results/${CLEAN_SPECIES}_${CLEAN_SEROTYPE}_assembled_results.tsv" > "EB_blast_results/filtered_${CLEAN_SPECIES}_${CLEAN_SEROTYPE}_assembled_results.tsv"
done
fi
done
# Ensure serotypes are loaded and process results for each serotype
mapfile -t LINES < "$SPECIES_SEROTYPE_FILE"
for LINE in "${LINES[@]}"; do
IFS=$'\t' read -r SPECIES SEROTYPE <<< "$LINE"
CLEAN_SEROTYPE=$(echo "$SEROTYPE" | sed 's/ /_/g; s/-/_/g; s/\.//g; s/,//g')
CLEAN_SPECIES=$(echo "$SPECIES" | sed 's/ /_/g')
if [[ "$SEROTYPE" != "NA" ]]; then
GENOME_DIR="EB_genomes/${CLEAN_SPECIES}_${CLEAN_SEROTYPE}"
TOTAL_DRAFT_COUNT=$(esearch -db assembly -query "$SPECIES AND $SEROTYPE AND (latest[filter] AND all[filter] NOT complete genome[filter])" | xtract -pattern Count -element Count)
else
CLEAN_SEROTYPE=""
GENOME_DIR="EB_genomes/${CLEAN_SPECIES}"
TOTAL_DRAFT_COUNT=$(esearch -db assembly -query "$SPECIES[Organism] AND (latest[filter] AND all[filter] NOT complete genome[filter])" | xtract -pattern Count -element Count)
fi
FILTERED_BLAST_OUTPUT="EB_blast_results/filtered_${CLEAN_SPECIES}_${CLEAN_SEROTYPE}_results.tsv"
TOTAL_COMPLETE_GENOMES=$(find "$GENOME_DIR" -type f -name "*_genomic.fna" | wc -l)
if [[ "$MODE" == "heavy" && "$TOTAL_COMPLETE_GENOMES" -lt 30 && "$SRA_FLAG" == "on" ]]; then
SRA_REQUESTED="$SRA_NUM"
TOTAL_ASSEMBLED_GENOMES=$(find "$ASSEMBLY_DIR/${CLEAN_SPECIES}_${CLEAN_SEROTYPE}" -mindepth 1 -type d | wc -l)
SRA_DROPPED=$((SRA_REQUESTED - TOTAL_ASSEMBLED_GENOMES))
else
SRA_REQUESTED=0
SRA_DROPPED=0
TOTAL_ASSEMBLED_GENOMES=0
fi
# Extract all unique gene IDs detected either in complete or assembled genomes
GENE_WITH_HITS=$(awk '{print $1}' "$FILTERED_BLAST_OUTPUT" 2>/dev/null; awk '{print $1}' "EB_blast_results/filtered_${CLEAN_SPECIES}_${CLEAN_SEROTYPE}_assembled_results.tsv" 2>/dev/null | sort -u)
mapfile -t ALL_GENES < <(grep "^>" "$GENE_FILE" | sed 's/>//' | awk '{print $1}')
for GENE_ID in "${ALL_GENES[@]}"; do
if grep -q "^${GENE_ID}$" <<< "$GENE_WITH_HITS"; then
COMPLETE_GENOMES_WITH_TARGET_GENES=$(awk -v gene="$GENE_ID" '$1 == gene {print $2}' "$FILTERED_BLAST_OUTPUT" | sort -u | wc -l)
if [[ "$MODE" == "heavy" && "$TOTAL_COMPLETE_GENOMES" -lt 30 && "$SRA_FLAG" == "on" ]]; then
ASSEMBLED_GENOMES_WITH_TARGET_GENES=$(awk -v gene="$GENE_ID" '$1 == gene {print $2}' "EB_blast_results/filtered_${CLEAN_SPECIES}_${CLEAN_SEROTYPE}_assembled_results.tsv" 2>/dev/null | awk -F'_' '{print $1"_"$2}' | sort -u | wc -l)
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
echo -e "$SPECIES,$SEROTYPE,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_COUNT,$TOTAL_COMPLETE_GENOMES,$SRA_REQUESTED,$SRA_DROPPED,$TOTAL_ASSEMBLED_GENOMES,$COMPLETE_GENOMES_WITH_TARGET_GENES,$ASSEMBLED_GENOMES_WITH_TARGET_GENES,$PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES%,$PERCENT_WITH_TARGET_GENES_ASSEMBLED_GENOMES%" >> "$OUTPUT_FILE"
else
echo -e "$SPECIES,$SEROTYPE,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_COUNT,$TOTAL_COMPLETE_GENOMES,$COMPLETE_GENOMES_WITH_TARGET_GENES,$PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES%"  >> "$OUTPUT_FILE"
fi
else
# If no hits were found for the gene in a serotype, output is recorded as 0%
if [[ "$MODE" == "heavy" ]]; then
echo -e "$SPECIES,$SEROTYPE,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_COUNT,$TOTAL_COMPLETE_GENOMES,$SRA_REQUESTED,$SRA_DROPPED,$TOTAL_ASSEMBLED_GENOMES,0,0,0%,0%" >> "$OUTPUT_FILE"
else
echo -e "$SPECIES,$SEROTYPE,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_COUNT,$TOTAL_COMPLETE_GENOMES,0,0%" >> "$OUTPUT_FILE"
fi
fi
done
done
echo "Analysis complete. Results saved in $OUTPUT_FILE"
