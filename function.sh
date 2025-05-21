#!/bin/bash
MONOPHASIC_TYPHIMURIUM_LIST="$(pwd)/database/monophasic_Typhimurium_list.txt"
ASSEMBLY_LEVEL="complete"
# Build a function to download genus-level genomes using ncbi-genome-download
function download_genus() {
local genus="$1"
local output_dir="$2"
local genus_dir="${output_dir}/${genus}/aggregated"
if [[ -n "$(find "$genus_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for $genus, skipping"
return
fi
mkdir -p "$genus_dir"
echo "Downloading genus: $genus"
ncbi-genome-download bacteria --genera "$genus" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$genus_dir"
find "$genus_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$0" && mv "${0%.gz}" "'"$genus_dir"'"' {} \;
rm -rf "$genus_dir/genbank"
echo "Downloaded and organized genomes for $genus"
}
# Build a function to get all species names under each target genus
function get_species_list() {
local target_genus="$1"
local output_dir="$2"
mkdir -p "$output_dir"
esearch -db taxonomy -query "${target_genus}[Subtree]" | efetch -format xml | xtract -pattern Taxon -element ScientificName | awk 'NF == 2 && $1 ~ /^[A-Z][a-z]+$/ && $2 ~ /^[a-z]+$/ {print $0}' | \
while read -r SPECIES; do
echo "$SPECIES" >> "$GENOME_DIR/species_list.txt"
done
echo "Saved all species names of $target_genus in species_list.txt"
}

# Build a function to download species-level genomes using ncbi-genome-download
function download_species() {
local target_genus="$1"
local output_dir="$2"
local species_list_file="${output_dir}/species_list.txt"
while read -r SPECIES; do
local genus=$(echo "$SPECIES" | awk '{print $1}')
if [[ "$genus" != "$target_genus" ]]; then
continue
fi
local clean_species=$(echo "$SPECIES" | sed 's/ /_/g')
local species_dir="${output_dir}/${genus}/${clean_species}"
[[ "$SPECIES" == "Salmonella enterica" ]] && species_dir="$species_dir/aggregated"
mkdir -p "$species_dir"
if [[ -n "$(find "$species_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for $SPECIES, skipping donwloading"
continue
fi
echo "Downloading genomes for $SPECIES"
ncbi-genome-download bacteria --genera "$SPECIES" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$species_dir" --verbose
# Move genomic FASTA files to correct location
if [[ -d "$species_dir/genbank" ]]; then
find "$species_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$0" && mv "${0%.gz}" "'"$species_dir"'"' {} \;
rm -rf "$species_dir/genbank"
fi
if [[ -z "$(find "$species_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "No genomes found for $SPECIES. Remove empty directory."
rm -rf "$species_dir"
fi
done < "$species_list_file"
echo "Downloaded and organized species genomes for $target_genus"
}
# Build a function to get all subspecies names under Salmonella enterica
function get_salmonella_subsp_list() {
local output_dir="$1"
mkdir -p "$output_dir"
esearch -db taxonomy -query "Salmonella enterica[Subtree]" | efetch -format xml | \
xtract -pattern Taxon -element ScientificName | grep -v -E "serovar|str\." | sort -u | grep -v "Salmonella enterica$" | sed 's/Salmonella enterica subsp. //' > "${output_dir}/salmonella_subspecies_list.txt"
echo "Saved all Salmonella enterica subsp. names in salmonella_subspecies_list.txt"
}

function download_salmonella_subsp() {
local output_dir="$1"
local subspecies_list="${output_dir}/salmonella_subspecies_list.txt"
echo "Downloading Salmonella enterica subspecies"
while read -r subspecies; do
local subspecies_dir="$GENOME_DIR/Salmonella/Salmonella_enterica/$subspecies"
[[ "$subspecies" == "enterica" ]] && subspecies_dir="$subspecies_dir/aggregated"
mkdir -p "$subspecies_dir"
if [[ -n "$(find "$subspecies_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for *Salmonella enterica* subsp. $subspecies"
continue
fi
echo "Downloading genomes for Salmonella enterica subsp. $subspecies"
ncbi-genome-download bacteria --genera "Salmonella enterica subsp. $subspecies" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$subspecies_dir"
if [[ -d "$subspecies_dir/genbank" ]]; then
find "$subspecies_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$0" && mv "${0%.gz}" "'"$subspecies_dir"'"' {} \;
rm -rf "$subspecies_dir/genbank"
fi
if [[ -z "$(find "$subspecies_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "No genomes found for $subspecies. Remove empty directory."
rm -rf "$subspecies_dir"
fi
done < "$subspecies_list"
}

# Build a function to get all serotype names under Salmonella enterica subsp. enterica
function get_salmonella_serotype_list() {
local output_dir="$1"
mkdir -p "$output_dir"
esearch -db taxonomy -query "Salmonella enterica subsp. enterica[Subtree]" | efetch -format xml | xtract -pattern Taxon -element ScientificName | grep -v -E "str\.|var\." | sort -u | grep -v "Salmonella enterica subsp. enterica$" > "${output_dir}/salmonella_serotype_list.txt"
echo "Saved all serotype names in salmonella_serotype_list.txt"
# Remove all Salmonella enterica subsp. enterica serovar prefix
sed -i 's/Salmonella enterica subsp. enterica serovar //g' "${output_dir}/salmonella_serotype_list.txt"
# Delete taxids that are already stored in all_monophasic_Typhimurium_taxids.txt
awk 'NR==FNR {monophasic[$0]; next} !($0 in monophasic)' "$MONOPHASIC_TYPHIMURIUM_LIST" "$GENOME_DIR/salmonella_serotype_taxids.txt" > tmp && mv tmp "${output_dir}/salmonella_serotype_taxids.txt"
} 

function download_salmonella_serotype() {
local output_dir="$1"
echo "Downloading all Salmonella enterica subsp. enterica serotype complete genomes..."
monophasic_dir="$GENOME_DIR/Salmonella/Salmonella_enterica/enterica/monophasic_Typhimurium"
mkdir -p "$monophasic_dir"
while read -r serotype; do
local serotype_dir="${GENOME_DIR}/Salmonella/Salmonella_enterica/enterica/${serotype}"
mkdir -p "$serotype_dir"
if [[ -n "$(find "$serotype_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for Salmonella enterica $serotype, skipping."
continue
fi
echo "Downloading genomes for Salmonella $serotype"
ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $serotype" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$serotype_dir" --verbose
if [[ -d "$serotype_dir/genbank" ]]; then
find "$serotype_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$1" && mv "${1%.gz}" "'"$serotype_dir"'"' _ {} \;
rm -rf "$serotype_dir/genbank"
fi
if [[ -z "$(find "$serotype_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "No genomes found for $subspecies. Remove empty directory."
rm -rf "$serotype_dir"
fi
done < "${output_dir}/salmonella_serotype_list.txt"
echo "Downloading all monophasic Typhimurium genomes"
while read -r monophasic; do
ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $monophasic" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$monophasic_dir" --verbose
if [[ -d "$monophasic_dir/genbank" ]]; then
find "$monophasic_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$1" && mv "${1%.gz}" "'"$monophasic_dir"'"' _ {} \;
rm -rf "$monophasic_dir/genbank"
fi
done < "$MONOPHASIC_TYPHIMURIUM_LIST"
}

function download_single_serotype() {
local serotype="$1"
local output_dir="${GENOME_DIR}/Salmonella/Salmonella_enterica/enterica/${serotype}"
echo "Downloading genomes for Salmonella $serotype..."
mkdir -p "$output_dir"
ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $serotype" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$output_dir" --verbose
if [[ -d "$output_dir/genbank" ]]; then
find "$output_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$1" && mv "${1%.gz}" "'"$output_dir"'"' _ {} \;
rm -rf "$output_dir/genbank"
fi
if [[ -z "$(find "$output_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "No genomes found for $subspecies. Remove empty directory."
rm -rf "$output_dir"
fi
}

function move_unclassified_genomes() {
local output_dir="$1"
find "$output_dir" -mindepth 1 -type d | while read -r level_dir; do
aggregated_dir="${level_dir}/aggregated"
unclassified_dir="${level_dir}/unclassified"
[[ ! -d "$aggregated_dir" ]] && continue
mkdir -p "$unclassified_dir"
declare -A classified_files=()
while IFS= read -r -d '' genome_path; do
genome_basename=$(basename "$genome_path")
classified_files["$genome_basename"]=1
done < <(find "$level_dir" -mindepth 1 -maxdepth 1 -type d \
! \( -name "aggregated" -o -name "unclassified" \) \
-exec find {} -type f -name "*_genomic.fna" -print0 \;)
while IFS= read -r -d '' genome_file; do
genome_basename=$(basename "$genome_file")
if [[ -z "${classified_files["$genome_basename"]}" ]]; then
echo "Moving unclassified genome: $genome_basename"
mv -f "$genome_file" "$unclassified_dir/" || echo "Error: Failed to move $genome_file"
fi
done < <(find "$aggregated_dir" -type f -name "*_genomic.fna" -print0)
rm -rf "$aggregated_dir"
done
echo "Organizing genomes into unclassified directories and removing aggregated directories"
}

function create_blastdb() {
local INPUT_FASTA="$1"
local DB_PATH="$2"
makeblastdb -in "$INPUT_FASTA" -dbtype nucl -out "${DB_PATH}"
echo "BLAST database created for ${DB_PATH}"
}

function build_blastdb() {
local input_dir="$1"
local output_dir="$2"
mkdir -p "$output_dir"
echo "Building BLAST database from genomes in $input_dir"
# Create BLAST database for all genus
find "$input_dir" -mindepth 1 -maxdepth 1 -type d | while read -r GENUS_DIR; do
GENUS=$(basename "$GENUS_DIR")
echo "Processing genus: $GENUS"
GENUS_FASTA="${GENUS_DIR}/${GENUS}_all_genomes.fna"
find "$GENUS_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$GENUS_FASTA"
create_blastdb "$GENUS_FASTA" "${output_dir}/${GENUS}"
# Species-level databases
find "$GENUS_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r SPECIES_DIR; do
SPECIES=$(basename "$SPECIES_DIR")
[[ "$SPECIES" == "unclassified" ]] && continue
SPECIES_FASTA="${SPECIES_DIR}/${SPECIES}_all_genomes.fna"
find "$SPECIES_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$SPECIES_FASTA"
create_blastdb "$SPECIES_FASTA" "${output_dir}/${SPECIES}"
if [[ "$SPECIES" == "Salmonella_enterica" ]]; then
find "$SPECIES_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r SUBSPECIES_DIR; do
SUBSPECIES=$(basename "$SUBSPECIES_DIR")
[[ "$SUBSPECIES" == "unclassified" ]] && continue
SUBSPECIES_FASTA="${SUBSPECIES_DIR}/Salmonella_enterica_subsp_${SUBSPECIES}_all_genomes.fna"
find "$SUBSPECIES_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$SUBSPECIES_FASTA"
create_blastdb "$SUBSPECIES_FASTA" "${output_dir}/Salmonella_enterica_subsp_${SUBSPECIES}"
if [[ "$SUBSPECIES" == "enterica" ]]; then
find "$SUBSPECIES_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r SEROTYPE_DIR; do
SEROTYPE=$(basename "$SEROTYPE_DIR")
[[ "$SEROTYPE" == "unclassified" ]] && continue
SEROTYPE_FASTA="${SEROTYPE_DIR}/${SEROTYPE}_genomes.fna"
find "$SEROTYPE_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$SEROTYPE_FASTA"
create_blastdb "$SEROTYPE_FASTA" "${output_dir}/Salmonella_${SEROTYPE}"
done
fi
done
fi
done
done
echo "Finished building BLAST databases"
}

function extract_taxon_info() {
local input="$1"
local taxon_name=""
read -ra words <<< "$input"
if [[ "${words[1]}" == "monophasic" ]]; then
taxon_name="Salmonella enterica subsp. enterica serovar monophasic Typhimurium"
elif [[ "${words[1]}" =~ ^[A-Z] || "${words[1]}" =~ ":" ]]; then
taxon_name="Salmonella enterica subsp. enterica serovar ${words[1]}"
else
taxon_name="$input"
fi
echo "$taxon_name"
}

function get_total_genomes_count() {
local input="$1"
local assembly_level="$2"
local total_genomes=0
local query="$(extract_taxon_info "$input")"
if [[ "$query" == "Salmonella enterica subsp. enterica serovar monophasic Typhimurium" ]]; then
while read -r actual_name; do
count=$(( $(ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $actual_name" --assembly-level "$assembly_level" --section genbank --dry-run  | wc -l) - 1 ))
((total_genomes += count))
done < "$MONOPHASIC_TYPHIMURIUM_LIST"
else
total_genomes=$(( $(ncbi-genome-download bacteria --genera "$query" --assembly-level "$assembly_level" --section genbank --dry-run --verbose | wc -l) - 1 ))
total_genomes=$(( total_genomes < 0 ? 0 : total_genomes ))
fi
echo "$total_genomes"
}

function calculate_sample_size_and_iterations(){
local total_genomes="$1"
local max_iterations=20
if [[ "$total_genomes" -eq 0 ]]; then
echo "0 0"
return
fi
if [[ "$total_genomes" -le 100 ]]; then
echo "$total_genomes 1"
return
fi
# Cochran's formula constants
local Z=1.96
local p=0.5
local e=0.05
local n0=$(echo "scale=6; ($Z^2 * $p * (1 - $p)) / ($e^2)" | bc -l)
# Finite population correction (FPC)
local sample_size=$(echo "scale=6; ($n0 * $total_genomes) / ($total_genomes + $n0 -1)" | bc -l)
sample_size=$(echo "$sample_size" | awk '{printf "%.0f", $1}')
# Compute number of iterations using square-root scaling
if [[ "$sample_size" -gt 0 ]]; then
local iterations=$(echo "scale=6; sqrt($total_genomes / (2 * $sample_size))" | bc -l)
iterations=$(echo "$iterations" | awk '{printf "%.0f", $1}')
else
local iterations=1
fi
if [[ "$iterations" -lt 1 ]]; then
iterations=1
elif [[ "$iterations" -gt "$max_iterations" ]]; then
iterations="$max_iterations"
fi
echo "$sample_size $iterations"
}

function download_random_draft_genomes() {
local input="$1"
local sample_size="$2"
local output_dir="$3"
local iteration="$4"
local query="$(extract_taxon_info "$input")"
local taxon_dir_name="${query// /_}"
taxon_dir_name="${taxon_dir_name//./}"
local iteration_dir="${output_dir}/${taxon_dir_name}/genomes_${iteration}"
mkdir -p "$iteration_dir"
local accessions=""
if [[ "$query" == "Salmonella enterica subsp. enterica serovar monophasic Typhimurium" ]]; then
while read -r actual_name; do
accessions+=$(ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $actual_name" --assembly-level contig --section genbank --dry-run | tail -n +2 | awk -F '/' '{print $NF}' | awk '{print $1}')
done < "$MONOPHASIC_TYPHIMURIUM_LIST"
else
accessions=$(ncbi-genome-download bacteria --genera "$query" --assembly-level contig --section genbank --dry-run | tail -n +2 | awk -F '/' '{print $NF}' | awk '{print $1}')
fi
valid_accessions=$(echo "$accessions" | grep -E '^GCA_[0-9]+\.[0-9]+$' | shuf -n "$sample_size")
echo "$valid_accessions" > "${iteration_dir}/selected_accessions.txt"
ncbi-genome-download bacteria --assembly-accessions "$iteration_dir/selected_accessions.txt" --formats fasta --assembly-level contig --section genbank  --output-folder "$iteration_dir"
find "$iteration_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$0" && mv "${0%.gz}" "'"$iteration_dir"'"' {} \;
rm -rf "$iteration_dir/genbank"
} 

function perform_blast(){
local query_gene="$1"
local perc_identity="$2"
local output_dir="$3"
local iteration="$4" # only used for draft genomes
local taxon="$5"
if [[ -n "$iteration" ]]; then
echo "Processing draft genomes for iteration $iteration: $taxon"
local standard_taxon="$(extract_taxon_info "$taxon")"
local safe_taxon="${standard_taxon// /_}"
safe_taxon="${safe_taxon//./}"
local genome_dir="${DRAFT_GENOMES_DIR}/${safe_taxon}/genomes_${iteration}"
local blast_db="${DRAFT_BLAST_DB_DIR}/${safe_taxon}/iteration_${iteration}"
mkdir -p "$(dirname "$blast_db")"
if [ ! -d "$genome_dir" ] || [ -z "$(find "$genome_dir" -type f -name "*_genomic.fna" 2>/dev/null)" ]; then
echo "Error: No genomic.fna files found in $genome_dir. Check download process." >&2
exit 1
fi
local concatenated_genome="${genome_dir}/combined.fna"
echo "Concatenating genome FASTAs for $taxon..."
find "$genome_dir" -name "*_genomic.fna" -exec cat {} + > "$concatenated_genome"
echo "Building BLAST DB for $taxon..."
makeblastdb -in "$concatenated_genome" -dbtype nucl -out "$blast_db"
local blast_output="${output_dir}/${taxon}/iteration_${iteration}_draft_blast_results.txt"
mkdir -p "$(dirname "$blast_output")"
echo "Running BLAST for $taxon..."
blastn -query "$query_gene" -db "$blast_db" -out "$blast_output" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
-perc_identity "$perc_identity" -max_target_seqs 7000
echo "BLAST results saved to: $blast_output"
else
# Skip empty lines
[[ -z "$taxon" ]] && return
local blast_db_name="${taxon// /_}"
blast_db_name="${blast_db_name//./}"
local blast_db="${BLAST_DB_DIR}/${blast_db_name}"
local clean_taxon="$(extract_taxon_info "$taxon")"
local blast_output_name="${clean_taxon// /_}"
blast_output_name="${blast_output_name//./}"
local blast_output="${output_dir}/${blast_output_name}_complete_blast_results.txt"
blastn -query "$query_gene" -db "$blast_db" -out "$blast_output" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
-perc_identity "$perc_identity" -max_target_seqs 7000
echo "BLAST results saved to $blast_output"
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
local filtered_result="${output_dir}/${taxon}/filtered_iteration_${iteration}_draft_blast_results.txt"
mkdir -p "$(dirname "$filtered_result")"
awk -v cov="$coverage_threshold" '($13 > 0) && (($4 / $13 * 100) >= cov) {print $0}' "$blast_result_file" > "$filtered_result"
echo "Filtered BLAST results for draft genomes are saved to $filtered_result"
else
local filtered_result="$output_dir/filtered_${base_name}"
awk -v cov="$coverage_threshold" '($13 > 0) && (($4 / $13 * 100) >= cov) {print $0}' "$blast_result_file" > "$filtered_result"
echo "Filtered BLAST results for complete genomes are saved to $filtered_result"
fi
}
