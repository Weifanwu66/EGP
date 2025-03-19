#!/bin/bash
# set initial variables and paths
MODE="light"
MIN_COVERAGE=80
MIN_IDENTITY=90
WORK_DIR=$(pwd)
TAXON_FILE=""
GENE_FILE=""
OUTPUT_FILE="${WORK_DIR}/result/gene_summary.csv"
EB_DRAFT_GENOMES_DIR="${WORK_DIR}/database/draft_genomes"
BLAST_DB_DIR="${WORK_DIR}/database/complete_blast_db"
BLAST_RESULT_DIR="${WORK_DIR}/result/complete_blast_results"
FILTERED_BLAST_RESULT_DIR="${WORK_DIR}/result/filtered_complete_blast_results"
DRAFT_BLAST_DB_DIR="${WORK_DIR}/database/draft_blast_db"
DRAFT_BLAST_RESULT_DIR="${WORK_DIR}/result/draft_blast_results"
FILTERED_DRAFT_BLAST_RESULT_DIR="${WORK_DIR}/result/filtered_draft_blast_results"
OVERWRITE=false
# usage setup
usage(){
echo "Usage: $0 -g GENE_FILE [-t TAXON] [-c COVERAGE] [-i IDENTITY] [-o OUTPUT_FILE] [--draft-sample N] [--mode light|heavy]"
echo "-g GENE_FILE : FASTA file with target gene sequence (required)." 
echo "-t TAXON_FILE : File containing target taxons (one target per line) or a single taxon name (e.g., "Salmonella Typhimurium"); A taxon_file is always recommended and must be provided in heavy mode."
echo "-c COVERAGE : The minimum genome coverage (default: 80%)." 
echo "-i IDENTITY : The minimum percentage of identity (default: 90%)."
echo "-o OUTPUT_FILE : Output result file (default: EB_gene_summary.tsv)."
echo "--mode light|heavy : Run in light mode (complete genomes only) or heavy mode (complete+draft genomes, default: light)."
echo "--overwrite : If overwrite mode is enabled, previous results will be cleared (default: false)."
}
# Parse argument
while [[ "$#" -gt 0 ]]; do
case "$1" in
-t) TAXON_FILE="$2"; shift ;;
-g) GENE_FILE="$2"; shift ;;
-c) MIN_COVERAGE="$2"; shift ;;
-i) MIN_IDENTITY="$2"; shift ;;
-o) OUTPUT_FILE="$2"; shift ;;
--mode) MODE="$2"; shift ;;
--overwrite) OVERWRITE=true ;;
-h|--help) usage; exit 0 ;;
*) echo "Invalid option: $1"; usage; exit 1 ;;
esac
shift
done 
if [[ "$OVERWRITE" == true ]]; then
rm -rf "$BLAST_RESULT_DIR" "$FILTERED_BLAST_RESULT_DIR" "$DRAFT_BLAST_RESULT_DIR" "$FILTERED_BLAST_RESULT_DIR" "$OUTPUT_FILE"
fi
# Ensure gene file is provided
if [[ -z "$GENE_FILE" ]]; then
echo "Error! Please provide a gene sequence file!"
usage
exit 1
fi
# Ensure if heavy mode is enabled, a taxon file must be provided
if [[ "$MODE" == "heavy" && -z "$TAXON_FILE" ]]; then
echo "Error: A taxon file is required in heavy mode."
exit 1
fi
# Output directory setup
mkdir -p "$BLAST_RESULT_DIR" "$FILTERED_BLAST_RESULT_DIR"
source "${WORK_DIR}/function.sh" || { echo "Error sourcing function.sh". exit 1; }
# Set delimiter as a space
if [[ -n "$TAXON_FILE" ]]; then
DELIMITER=" "
sed 's/[\t,]\+/ /g' "$TAXON_FILE" > "${TAXON_FILE}_processed"
# Ensure if a taxon file is provided, the genus should be one of the target Enterobacteriaceae
awk -v delim="$DELIMITER" '
BEGIN { FS=delim; OFS=delim }
{ 
gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
genus = $1
species_or_serotype = ($2 != "") ? $2 : ""
if (genus !~ /^(Proteus|Salmonella|Escherichia|Citrobacter|Enterobacter|Klebsiella|Shigella|Cronobacter)$/) {
print "Error: Invalid genus in taxon file. Allowed: Proteus, Salmonella, Escherichia, Citrobacter, Enterobacter, Klebsiella, Shigella, Cronobacter."
print "Your line:", $0
exit 1 
}
}' "${TAXON_FILE}_processed" || exit 1
rm "${TAXON_FILE}_processed"
fi
# start with core processing functions
process_complete_genomes() {
local taxon="$1"
echo "Processing complete genomes for $taxon"
local blast_db_name="${taxon// /_}"
blast_db_name="${blast_db_name//./}"
local blast_output="${BLAST_RESULT_DIR}/${blast_db_name}_complete_blast_results.txt"
if [[ ! -s "$blast_output" ]]; then
perform_blast "$GENE_FILE" "$MIN_IDENTITY" "$BLAST_RESULT_DIR" "" "$taxon"
echo "BLAST result saved to $blast_output"
else
echo "Skipping BLAST: results already exist for $taxon"
fi
for blast_result_file in "$blast_output"; do
filter_blast_results "$blast_result_file" "$FILTERED_BLAST_RESULT_DIR" "$MIN_COVERAGE" "complete"
done
echo "Finished processing complete genomes for $taxon"
}

process_draft_genomes() {
local taxon="$1"
local total_draft_genomes=$(get_total_genomes_count "$taxon" "contig")
read -r sample_size iterations <<< "$(calculate_sample_size_and_iterations "$total_draft_genomes")"
echo "Processing $taxon | Total draft genomes: $total_draft_genomes. Running $iterations iterations (max 20)."
for ((i=1; i<=iterations; i++)); do
echo "Starting iterations $i/$iterations for $taxon"
download_random_draft_genomes "$taxon" "$sample_size" "$EB_DRAFT_GENOMES_DIR" "$i"
perform_blast "$GENE_FILE" "$MIN_IDENTITY" "$DRAFT_BLAST_RESULT_DIR" "$i" ""
local blast_db_name="${taxon// /_}"
blast_db_name="${blast_db_name//./}"
local blast_result_file="${DRAFT_BLAST_RESULT_DIR}/${blast_db_name}/iteration${i}_draft_blast_results.txt"
filter_blast_results "$blast_result_file" "$FILTERED_DRAFT_BLAST_RESULT_DIR" "$MIN_COVERAGE" "draft"
done
}
# Set up the output file headers
if [[ "$MODE" == "heavy" ]]; then
echo -e "Organism,Gene_ID,Min_percentage_of_coverage,Min_percentage_of_identity,Total_draft_genomes,Total_complete_genomes,Draft_genomes_sample_size,Number_of_iterations,Complete_genomes_with_target_genes,Draft_genomes_with_target_genes,Percentage_with_target_genes_complete_genomes,Percentage_with_target_genes_draft_genomes" > "$OUTPUT_FILE"
else
echo -e "Organism,Gene_ID,Min_percentage_of_coverage,Min_percentage_of_identity,Total_draft_genomes,Total_complete_genomes,Complete_genomes_with_target_genes,Percentage_with_target_genes_complete_genomes" > "$OUTPUT_FILE"
fi
# Initialize TAXON_LIST
TAXON_LIST=""
if [[ -n "$TAXON_FILE" ]]; then
if [[ ! -f "$TAXON_FILE" ]]; then
tmp_taxon_file=$(mktemp)
echo "$TAXON_FILE" > "$tmp_taxon_file"
TAXON_FILE="$tmp_taxon_file"
trap "rm -f '$tmp_taxon_file'" EXIT
fi
TAXON_LIST=$(< "$TAXON_FILE")
else
# Rmove numeric suffix and convert underscores and dots to spaces
TAXON_LIST=$(find "$BLAST_DB_DIR" -name "*.nsq" -exec basename {} .nsq \; | sed -E '/[._][0-9]{2}$/s/[._][0-9]{2}$//' | sort -u |  sed 's/_/ /g' | sed -E 's/subsp /subsp. /g')
fi
# Process each taxon in TAXON_LIST
while IFS= read -r taxon; do
[[ -z "$taxon" ]] && continue
echo "Processing $taxon in $MODE mode"
process_complete_genomes "$taxon"
[[ "$MODE" == "heavy" ]] && process_draft_genomes "$taxon"
TOTAL_COMPLETE_GENOMES=$(get_total_genomes_count "$taxon" "complete") 
TOTAL_DRAFT_GENOMES=$(get_total_genomes_count "$taxon" "contig")
local_taxon="${taxon// /_}"
local_taxon="${local_taxon//./}"
if [[ "$MODE" == "heavy" && "$TOTAL_DRAFT_GENOMES" -gt 0 ]]; then
read -r DRAFT_SAMPLE_SIZE ITERATIONS <<< "$(calculate_sample_size_and_iterations "$TOTAL_DRAFT_GENOMES")"
else
DRAFT_SAMPLE_SIZE=0
ITERATIONS=0
fi
GENE_WITH_HITS=$(awk '{print $1}' "$FILTERED_BLAST_RESULT_DIR/filtered_${local_taxon}_complete_blast_results.txt" 2>/dev/null)
[[ "$MODE" == "heavy" ]] && GENE_WITH_HITS+=$'\n'$(awk '{print $1}' "$FILTERED_DRAFT_BLAST_RESULT_DIR/$local_taxon"/* 2>/dev/null)
GENE_WITH_HITS=$(echo "$GENE_WITH_HITS" | sort -u)
# Process all genes in query gene file
mapfile -t ALL_GENES < <(grep "^>" "$GENE_FILE" | sed 's/>//' | awk '{print $1}')
for GENE_ID in "${ALL_GENES[@]}"; do
if grep -qx "$GENE_ID" <<< "$GENE_WITH_HITS"; then
COMPLETE_GENOMES_WITH_TARGET_GENES=$(awk -v gene="$GENE_ID" '$1 == gene {print $2}' "$FILTERED_BLAST_RESULT_DIR/filtered_${local_taxon}_complete_blast_results.txt" 2>/dev/null | sort -u | wc -l)
if [[ "$MODE" == "heavy" && "$TOTAL_DRAFT_GENOMES" -gt 0 ]]; then
TOTAL_DRAFT_GENOMES_WITH_TARGET_GENES=0
for ((i=1; i<="$ITERATIONS"; i++)); do
DRAFT_GENOMES_WITH_TARGET_GENES=$(awk -v gene="$GENE_ID" '$1 == gene {print $2}' "$FILTERED_DRAFT_BLAST_RESULT_DIR/$local_taxon/filtered_iteration${i}_draft_blast_results.txt" 2>/dev/null | awk -F'_' '{print substr($2, 1, 8)}' | sort -u | wc -l)
TOTAL_DRAFT_GENOMES_WITH_TARGET_GENES=$((TOTAL_DRAFT_GENOMES_WITH_TARGET_GENES + DRAFT_GENOMES_WITH_TARGET_GENES))
done
AVERAGE_DRAFT_GENOMES_WITH_TARGET_GENES=$(echo "scale=2; $TOTAL_DRAFT_GENOMES_WITH_TARGET_GENES/$ITERATIONS" | bc 2>/dev/null)
PERCENT_WITH_TARGET_GENES_DRAFT_GENOMES=$(echo "scale=2; ($AVERAGE_DRAFT_GENOMES_WITH_TARGET_GENES/$DRAFT_SAMPLE_SIZE) * 100" | bc 2>/dev/null)
else
DRAFT_GENOMES_WITH_TARGET_GENES=0
PERCENT_WITH_TARGET_GENES_DRAFT_GENOMES=0
fi

if [[ "$TOTAL_COMPLETE_GENOMES" -gt 0 ]]; then
PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES=$(echo "scale=2; ($COMPLETE_GENOMES_WITH_TARGET_GENES/$TOTAL_COMPLETE_GENOMES) * 100" | bc 2>/dev/null)
else
PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES=0
fi
if [[ "$MODE" == "heavy" ]]; then
echo -e "$taxon,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$DRAFT_SAMPLE_SIZE,$ITERATIONS,$COMPLETE_GENOMES_WITH_TARGET_GENES,$DRAFT_GENOMES_WITH_TARGET_GENES,${PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES}%,${PERCENT_WITH_TARGET_GENES_DRAFT_GENOMES}%" >> "$OUTPUT_FILE"
else
echo -e "$taxon,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$COMPLETE_GENOMES_WITH_TARGET_GENES,${PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES}%"  >> "$OUTPUT_FILE"
fi
else
# If no hits were found for the gene in a serotype, output is recorded as 0%
if [[ "$MODE" == "heavy" ]]; then
echo -e "$taxon,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$DRAFT_SAMPLE_SIZE,$ITERATIONS,0,0,0%,0%" >> "$OUTPUT_FILE"
else
echo -e "$taxon,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,0,0%" >> "$OUTPUT_FILE"
fi
fi
done
done <<< "$TAXON_LIST"
echo "Analysis complete. Results saved in $OUTPUT_FILE"
