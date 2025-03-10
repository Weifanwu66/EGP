#!/bin/bash
# set initial variables and paths
MODE="light"
MIN_COVERAGE=80
MIN_IDENTITY=90
WORK_DIR=$(pwd)
TAXON_FILE=""
GENE_FILE=""
SEROTYPE_TAXID_FILE="${WORK_DIR}/database/salmonella_serotype_taxids.txt"
OUTPUT_FILE="${WORK_DIR}/result/EB_gene_summary.csv"
EB_GENOMES_DIR="${WORK_DIR}/database/EB_complete_genomes"
EB_DRAFT_GENOMES_DIR="${WORK_DIR}/database/EB_draft_genomes"
BLAST_DB_DIR="${WORK_DIR}/database/EB_blast_db"
BLAST_RESULT_DIR="${WORK_DIR}/result/EB_blast_results"
FILTERED_BLAST_RESULT_DIR="${WORK_DIR}/result/filtered_EB_blast_results"
DRAFT_BLAST_DB_DIR="${WORK_DIR}/database/draft_EB_blast_db"
DRAFT_BLAST_RESULT_DIR="${WORK_DIR}/result/draft_EB_blast_results"
FILTERED_DRAFT_BLAST_RESULT_DIR="${WORK_DIR}/result/filtered_draft_EB_blast_results"
DRAFT_SAMPLE_SIZE=100
MAX_ITERATIONS=50
# usage setup
usage(){
echo "Usage: $0 -g GENE_FILE [-t TAXON_FILE] [-c COVERAGE] [-i IDENTITY] [-o OUTPUT_FILE] [--draft-sample N] [--mode light|heavy]"
echo "-g GENE_FILE : FASTA file with target gene sequence (required)." 
echo "-t TAXON_FILE : File containing your target taxons (one target per line); Must be provided in heavy mode."
echo "-c COVERAGE : The minimum genome coverage (default: 80%)." 
echo "-i IDENTITY : The minimum percentage of identity (default: 90%)."
echo "-o OUTPUT_FILE : Output result file (default: EB_gene_summary.tsv)."
echo "--draft-sample N : Number of draft genomes to randomly select per iteration (default: 100)."
echo "--mode light|heavy : Run in light mode (complete genomes only) or heavy mode (complete+draft genomes, default: light)."
}
# Parse argument
while [[ "$#" -gt 0 ]]; do
case "$1" in
-t) TAXON_FILE="$2"; shift ;;
-g) GENE_FILE="$2"; shift ;;
-c) MIN_COVERAGE="$2"; shift ;;
-i) MIN_IDENTITY="$2"; shift ;;
-o) OUTPUT_FILE="$2"; shift ;;
--draft-sample) DRAFT_SAMPLE_SIZE="$2"; shift ;;
--mode) MODE="$2"; shift ;;
-h|--help) usage; exit 0 ;;
*) echo "Invalid option: $1"; usage; exit 1 ;;
esac
shift
done 
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
source "${WORK_DIR}/function.sh"
# Auto detect delimiter of user's input file
if [[ -n "$TAXON_FILE" ]]; then
FIRST_LINE=$(head -n 1 "$TAXON_FILE")
if [[ "$FIRST_LINE" == *$'\t'* ]]; then
DELIMITER=$'\t'
elif [[ "$FIRST_LINE" == *" "* ]]; then
DELIMITER=$' '
elif [[ "$FIRST_LINE" == *","* ]]; then
DELIMITER=","
else
echo "Unknown delimiter in file. Ensure it is tab, comma, or space-separate."
exit 1
fi
# Ensure if a taxon file is provided, the genus should be one of the target Enterobacteriaceae
awk -v delim="$DELIMITER" '
BEGIN { FS=delim }
{ 
if ($1 !~ /^(Salmonella|Escherichia|Citrobacter|Enterobacter|Klebsiella|Shigella|Cronobacter)$/) {
print "Error: Invalid genus in taxon file. Allowed: Salmonella, Escherichia, Citrobacter, Enterobacter, Klebsiella, Shigella, Cronobacter."
print "Your line:", $0
exit 1 
}
}' "$TAXON_FILE" || exit 1

# start with core processing functions
process_complete_genomes() {
perform_blast "$GENE_FILE" "$MIN_IDENTITY" "$BLAST_RESULT_DIR" "" "$TAXON_FILE"
for blast_result_file in "${BLAST_RESULT_DIR}"/*_complete_blast_results.txt; do
filter_blast_results "$blast_result_file" "$FILTERED_BLAST_RESULT_DIR" "$MIN_COVERAGE" "complete"
done
echo "Finished processing complete genomes"
}

process_draft_genomes() {
local genus="$1"
local species_or_serotype="$2"
echo "Processing draft genomes"
local total_draft_genomes=$(get_total_genomes_count "$genus" "$species_or_serotype" "contig")
local iterations=$(( total_draft_genomes / DRAFT_SAMPLE_SIZE ))
(( iterations < 1 )) && iterations=1
(( iterations > MAX_ITERATIONS )) && iterations=$MAX_ITERATIONS
local taxon="${genus}"
[[ -n "$species_or_serotype" ]] && taxon+="_${species_or_serotype}"
echo "Processing $taxon | Total draft genomes: $total_draft_genomes. Running $iterations iterations (max 100) with sample size $DRAFT_SAMPLE_SIZE."
for ((i=1; i<=iterations; i++)); do
echo "Starting iterations $i/$iterations for $taxon"
download_random_draft_genomes "$genus" "$species_or_serotype" "$DRAFT_SAMPLE_SIZE" "$EB_DRAFT_GENOMES_DIR" "$i"
perform_blast "$GENE_FILE" "$MIN_IDENTITY" "$DRAFT_BLAST_RESULT_DIR" "$i" ""
local blast_result_file="${DRAFT_BLAST_RESULT_DIR}/${taxon}/iteration${i}_draft_blast_results.txt"
filter_blast_results "${blast_result_file}" "$FILTERED_DRAFT_BLAST_RESULT_DIR" "$MIN_COVERAGE" "draft"
done
}
echo "Running in $MODE mode"
echo "Processing complete genomes"
process_complete_genomes
# When heavy mode is on
if [[ "$MODE" == "heavy" ]]; then
echo "Processing draft genomes"
while IFS="$DELIMITER" read -r genus species_or_serotype || [[ -n "$genus" ]]; do
[[ -z "$genus" ]] && continue
process_draft_genomes "$genus" "$species_or_serotype"
done < "$TAXON_FILE"
fi

# Set up the output file headers
if [[ "$MODE" == "heavy" ]]; then
echo -e "Genus,Species_or_Serotype,Gene_ID,Min_percentage_of_coverage,Min_percentage_of_identity,Total_draft_genomes,Total_complete_genomes,Draft_genomes_sample_size,Number_of_iterations,Complete_genomes_with_target_genes,Draft_genomes_with_target_genes,Percentage_with_target_genes_complete_genomes,Percentage_with_target_genes_draft_genomes" > "$OUTPUT_FILE"
else
echo -e "Genus,Species_or_Serotype,Gene_ID,Min_percentage_of_coverage,Min_percentage_of_identity,Total_draft_genomes,Total_complete_genomes,Complete_genomes_with_target_genes,Percentage_with_target_genes_complete_genomes" > "$OUTPUT_FILE"
fi
# Ensure serotypes are loaded and process results for each serotype
while IFS="$DELIMITER" read -r genus species_or_serotype || [[ -n "$genus" ]]; do
# skip empty lines
[[ -z "$genus" ]] && continue
TOTAL_COMPLETE_GENOMES=$(get_total_genomes_count "$genus" "$species_or_serotype" "complete") 
TOTAL_DRAFT_GENOMES=$(get_total_genomes_count "$genus" "$species_or_serotype" "contig")
# Extract all unique gene IDs detected either in complete or in draft genomes
local_taxon="$genus"
[[ -n "$species_or_serotype" ]] && local_taxon+="_${species_or_serotype}"
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
echo -e "$genus,$species_or_serotype,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$DRAFT_SAMPLE_SIEZ,$ITERATIONS,$COMPLETE_GENOMES_WITH_TARGET_GENES,$DRAFT_GENOMES_WITH_TARGET_GENES,${PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES}%,${PERCENT_WITH_TARGET_GENES_DRAFT_GENOMES}%" >> "$OUTPUT_FILE"
else
echo -e "$genus,$species_or_serotype,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$COMPLETE_GENOMES_WITH_TARGET_GENES,${PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES}%"  >> "$OUTPUT_FILE"
fi
else
# If no hits were found for the gene in a serotype, output is recorded as 0%
if [[ "$MODE" == "heavy" ]]; then
echo -e "$genus,$species_or_serotype,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$DRAFT_SAMPLE_SIZE,$ITERATIONS,0,0,0%,0%" >> "$OUTPUT_FILE"
else
echo -e "$genus,$species_or_serotype,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,0,0%" >> "$OUTPUT_FILE"
fi
fi
done
done < "$TAXON_FILE"
echo "Analysis complete. Results saved in $OUTPUT_FILE"
