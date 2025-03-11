#!/bin/bash
function clean_taxons() {
echo "$1" | sed 's| |_|g; s|\.|_|g; s|-|_|g; s|,||g; s|\[||g; s|\]||g; s|/|_|g'
}

function extract_taxon_levels_from_directory() {
local name="$1"
local genus=""
local species="NA"
local serotype="NA"
IFS="_" read -r -a parts <<< "$name"
genus="${parts[0]}"
local second_part="${name#*_}"
if [[ -n "$second_part" ]]; then
if [[ "$second_part" =~ [A-Z] || "$second_part" =~ ":" ]]; then
species="enterica"
serotype="$second_part"
else
species="$second_part"
fi
fi
echo "$genus" "$species" "$serotype"
}

function  get_salmonella_serotype_taxid() {
local serotype="$1"
local taxid=$(awk -v s="$serotype" -F'\t' '$2 == s {print $1}' "$(pwd)/database/salmonella_serotype_taxids.txt")
echo "$taxid"
}

function get_total_genomes_count() {
local genus="$1"
local species_or_serotype="$2"
local assembly_level="$3"
if [[ "$species_or_serotype" =~ [A-Z] || "$species_or_serotype" =~ : ]]; then
taxid=$(get_salmonella_serotype_taxid "$species_or_serotype")
for TAXID in $taxid; do
COUNT=$(datasets summary genome taxon "$TAXID" --assembly-level "$assembly_level" 2>/dev/null | jq -r '.total_count')
total_genomes=$((total_genomes + COUNT))
done
else
if [[ -z "$species_or_serotype" ]]; then
taxon="$genus"
else
taxon="$genus $species_or_serotype"
fi
total_genomes=$(datasets summary genome taxon "$taxon" --assembly-level "$assembly_level" 2>/dev/null | jq -r '.total_count')
fi
echo "$total_genomes"
}

function download_random_draft_genomes() {
local genus="$1"
local species_or_serotype="$2"
local sample_size="$3"
local output_dir="$4"
local iteration="$5"
local taxon
if [[ -z "$species_or_serotype" ]]; then
taxon="$genus"
else
taxon="${genus}_${species_or_serotype}"
fi
local iteration_dir="${output_dir}/${taxon}/genomes_${iteration}"
mkdir -p "$iteration_dir"
local total_genomes
total_genomes=$(get_total_genomes_count "$genus $species_or_serotype" "draft")
[[ "$sample_size" -gt "$total_genomes" ]] && sample_size="$total_genomes"
local accessions
if [[ "$species_or_serotype" =~ [A-Z] || "$species_or_serotype" =~ : ]]; then
local taxids=$(get_salmonella_serotype_taxid "$species_or_serotype")
accessions=$(for TAXID in $taxids; do datasets summary genome taxon "$TAXID" --assembly-level contig | jq -r '.reports[].accession' done | shuf -n "$sample_size")
else
local taxon_query="$genus${species_or_serotype:+ $species_or_serotype}"
accessions=$(datasets summary genome taxon "$taxon_query" --assembly-level contig | jq -r '.reports[].accession' | shuf -n "$sample_size")
fi
echo "$accessions" > "${iteration_dir}/selected_accessions.txt"
datasets download genome accession --inputfile "${iteration_dir}/selected_accessions.txt" --filename "${iteration_dir}/download.zip"
unzip -q -o "${iteration_dir}/download.zip" -d "${iteration_dir}"
find "${iteration_dir}/ncbi_dataset" -type f -name "*_genomic.fna" -exec mv {} "${iteration_dir}" \;
rm "${iteration_dir}/download.zip"
rm -rf "${iteration_dir}/ncbi_dataset"
} 

function perform_blast(){
local query_gene="$1"
local perc_identity="$2"
local output_dir="$3"
local iteration="$4" # only used for draft genomes
local taxon_file="$5"  # only provide for complete genomes, always empty for draft genomes
if [[ -n "$iteration" ]]; then
echo "Processing draft genomes for iteration $iteration"
local taxon_dirs=("${EB_DRAFT_GENOMES_DIR}"/*/)
local taxon=$(basename "${taxon_dirs[0]%/}")
local genome_dir="${EB_DRAFT_GENOMES_DIR}/${taxon}/genomes_${iteration}"
local blast_db="${DRAFT_BLAST_DB_DIR}/${taxon}/iteration_${iteration}"
mkdir -p "$(dirname "$blast_db")"
local concatenated_genome="${genome_dir}/combined.fna"
cat "$genome_dir"/*_genomic.fna > "$concatenated_genome"
makeblastdb -in "$concatenated_genome" -dbtype nucl -out "$blast_db"
local blast_output="${output_dir}/${taxon}/iteration_${iteration}_draft_blast_results.txt"
mkdir -p "$(dirname "$blast_output")"
blastn -query "$query_gene" -db "$blast_db" -out "$blast_output" -outfmt 6 -perc_identity "$perc_identity"
echo "BLAST results saved to: $blast_output"
else
echo "Processing complete genomes"
if [[ -n "$taxon_file" ]]; then
while IFS="$DELIMITER" read -r genus species_or_serotype || [[ -n "$genus" ]]; do
# Skip empty lines
[[ -z "$genus" ]] && continue
local blast_db_name="${genus}"
[[ -n "$species_or_serotype" ]] && blast_db_name+="_${species_or_serotype// /_}"
local blast_db="${BLAST_DB_DIR}/${blast_db_name}"
local blast_output="${output_dir}/${blast_db_name}_complete_blast_results.txt"
blastn -query "$query_gene" -db "$blast_db" -out "$blast_output" -outfmt 6 -perc_identity "$perc_identity"
echo "BLAST results saved to $blast_output"
done < "$taxon_file"
else
echo "No taxon file provided. Processing all pre-built BLAST database"
find "$BLAST_DB_DIR" -name "*.nsq" -exec basename {} .nsq \; | sed -E '/[._][0-9]{2}$/s/[._][0-9]{2}$//' | sort -u | while read -r line; do
genus=$(echo "$line" | awk '{print $1}' | xargs)
species_or_serotype=$(echo "$line" | awk '{$1=""; print $0}' | xargs)
[[ "$genus" == "$line" ]] && species_or_serotype=""
local blast_db_name="$genus"
[[ -n "$species_or_serotype" ]] && blast_db_name+="_${species_or_serotype// /_}"
local blast_db="${BLAST_DB_DIR}/${blast_db_name}"
local blast_output="${output_dir}/${blast_db_name}_complete_blast_results.txt"
blastn -query "$query_gene" -db "$blast_db" -out "$blast_output" -outfmt 6 -perc_identity "$perc_identity"
echo "BLAST results saved to $blast_output"
done
fi
fi
echo "BLAST analysis completed. Results saved in: $output_dir"
}

function filter_blast_results() {
local blast_result_file="$1"  # Full path containing blast results
local output_dir="$2"
local coverage_threshold="$3"
local genome_type="$4"
local base_name="$(basename "$blast_result_file")"
if [[ "$genome_type" == "draft" ]]; then
local taxon="$(basename "$(dirname "$blast_result_file")")"
local iteration="$(grep -o '[0-9]\+' <<< "$base_name")"
local filtered_result="${output_dir}/${taxon}/filtered_iteration${iteration}_draft_blast_results.txt"
mkdir -p "$(dirname "$filtered_output")"
awk -v cov="$coverage_threshold" '(($4/($8 - $7 + 1)*100) >= cov) {print $0}' "$blast_result_file" > "$filtered_result"
echo "Filtered BLAST results for draft genomes are saved to $filtered_result"
else
local filtered_result="$output_dir/filtered_${base_name}"
awk -v cov="$coverage_threshold" '(($4/($8 - $7 + 1)*100) >= cov) {print $0}' "$blast_result_file" > "$filtered_result"
echo "Filtered BLAST results for complete genomes are saved to $filtered_result"
fi
}
