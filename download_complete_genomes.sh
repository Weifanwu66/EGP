#!/bin/bash
# Initial variables
WORKDIR=$(pwd)
GENOME_DIR="$WORKDIR/database/complete_genomes"
ASSEMBLY_LEVEL="complete"
MONOPHASIC_TYPHIMURIUM_LIST="$GENOME_DIR/monophasic_Typhimurium_list.txt"
# List of Enterobacteriaceae Genera
GENUS=("Proteus" "Escherichia" "Salmonella" "Shigella" "Klebsiella" "Enterobacter" "Cronobacter" "Citrobacter")
# Build a function to download genus-level genomes using ncbi-genome-download
download_genus() {
local genus="$1"
local genus_dir="$GENOME_DIR/$genus/aggregated"
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
get_species_list() {
local target_genus="$1"
esearch -db taxonomy -query "${target_genus}[Subtree]" | efetch -format xml | xtract -pattern Taxon -element ScientificName | awk 'NF == 2 && $1 ~ /^[A-Z][a-z]+$/ && $2 ~ /^[a-z]+$/ {print $0}' | \
while read -r SPECIES; do
echo "$SPECIES" >> "$GENOME_DIR/species_list.txt"
done
echo "Saved all species names of $target_genus in species_list.txt"
}
# Build a function to download species-level genomes using ncbi-genome-download
download_species() {
local target_genus="$1"
while read -r SPECIES; do
local genus=$(echo "$SPECIES" | awk '{print $1}')
if [[ "$genus" != "$target_genus" ]]; then
continue
fi
local clean_species=$(echo "$SPECIES" | sed 's/ /_/g')
local species_dir="$GENOME_DIR/$genus/$clean_species"
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
done < "$GENOME_DIR/species_list.txt"
echo "Downloaded and organized species genomes for $target_genus"
}
# Build a function to get all subspecies names under Salmonella enterica
get_salmonella_subsp_list() {
esearch -db taxonomy -query "Salmonella enterica[Subtree]" | efetch -format xml | \
xtract -pattern Taxon -element ScientificName | grep -v -E "serovar|str\." | sort -u | grep -v "Salmonella enterica$" | sed 's/Salmonella enterica subsp. //' > "$GENOME_DIR/salmonella_subspecies_list.txt"
echo "Saved all Salmonella enterica subsp. names in salmonella_subspecies_list.txt"
}
download_salmonella_subsp() {
echo "Downloading *Salmonella enterica* subspecies"
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
done < "$GENOME_DIR/salmonella_subspecies_list.txt"
}
# Build a function to get all serotype names under Salmonella enterica subsp. enterica
get_salmonella_serotype_list() {
esearch -db taxonomy -query "Salmonella enterica subsp. enterica[Subtree]" | efetch -format xml | xtract -pattern Taxon -element ScientificName | grep -v -E "str\.|var\." | sort -u | grep -v "Salmonella enterica subsp. enterica$" > "$GENOME_DIR/salmonella_serotype_list.txt"
echo "Saved all serotype names in salmonella_serotype_list.txt"
# Remove all Salmonella enterica subsp. enterica serovar prefix
sed -i 's/Salmonella enterica subsp. enterica serovar //g' "$GENOME_DIR/salmonella_serotype_list.txt"
# Delete taxids that are already stored in all_monophasic_Typhimurium_taxids.txt
awk 'NR==FNR {monophasic[$0]; next} !($0 in monophasic)' "$MONOPHASIC_TYPHIMURIUM_LIST" "$GENOME_DIR/salmonella_serotype_taxids.txt" > tmp && mv tmp "$GENOME_DIR/salmonella_serotype_taxids.txt"
} 
download_salmonella_serotype() {
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
done < "$GENOME_DIR/salmonella_serotype_list.txt"
echo "Downloading all monophasic Typhimurium genomes"
while read -r monophasic; do
ncbi-genome-download bacteria --genera "Salmonella enterica subsp. enterica serovar $monophasic" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$monophasic_dir" --verbose
if [[ -d "$monophasic_dir/genbank" ]]; then
find "$monophasic_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$1" && mv "${1%.gz}" "'"$monophasic_dir"'"' _ {} \;
rm -rf "$monophasic_dir/genbank"
fi
done < "$GENOME_DIR/monophasic_Typhimurium_list.txt"
}
# Function to categorize any genomes stored in aggregated directory but not present in their sibling directories as unclassified
move_unclassified_genomes() {
find "$GENOME_DIR" -mindepth 1 -type d | while read -r level_dir; do
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

# Download genomes 
for GENUS in "${GENUS[@]}"; do
download_genus "$GENUS"
get_species_list "$GENUS"
download_species "$GENUS" 
done
get_salmonella_subsp_list
download_salmonella_subsp
get_salmonella_serotype_list
download_salmonella_serotype
# Organize unclassified genomes
move_unclassified_genomes
