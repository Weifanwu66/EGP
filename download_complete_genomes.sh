#!/bin/bash
# Initial variables
WORKDIR=$(pwd)
GENOME_DIR="$WORKDIR/database/complete_genomes"
ASSEMBLY_LEVEL="complete"
MONOPHASIC_TYPHIMURIUM_TAXIDS="$GENOME_DIR/all_monophasic_Typhimurium_taxids.txt"
# List of Enterobacteriaceae Genera
EB_GENUS=("Escherichia" "Salmonella" "Shigella" "Klebsiella" "Enterobacter" "Cronobacter" "Citrobacter")
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
# Build a function to get all species taxid
get_species_taxids() {
for GENUS in "${EB_GENUS[@]}"; do
esearch -db taxonomy -query "${GENUS}[Subtree]" | efetch -format xml | xtract -pattern Taxon -element TaxId,ScientificName | awk 'NF == 3 && $2 ~ /^[A-Z][a-z]+$/ && $3 ~ /^[a-z]+$/ {print $1, $2, $3}' >> "$GENOME_DIR/species_taxids.txt"
done
echo "Saved all species Taxon IDs in species_taxids.txt"
}
# Build a function to download species-level genomes using ncbi-genome-download
download_species_by_taxid() {
while read -r TAXID SPECIES; do
local genus=$(echo "$SPECIES" | awk '{print $1}')
local clean_species=$(echo "$SPECIES" | sed 's/ /_/g')
local species_dir="$GENOME_DIR/$genus/$clean_species"
[[ "$SPECIES" == "Salmonella enterica" ]] && species_dir="$species_dir/aggregated"
mkdir -p "$species_dir"
if [[ -n "$(find "$species_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for $SPECIES, skipping donwloading"
continue
fi
echo "Downloading genomes for $SPECIES"
ncbi-genome-download bacteria --taxid "$TAXID" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$species_dir"
# Move genomic FASTA files to correct location
if [[ -d "$species_dir/genbank" ]]; then
find "$species_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$0" && mv "${0%.gz}" "'"$species_dir"'"' {} \;
rm -rf "$species_dir/genbank"
fi
if [[ -z "$(find "$species_dir" -maxdepth 1 -type d -name "genbank" 2>/dev/null)" ]]; then
echo "No genomes found for $SPECIES. Remove empty directory."
rm -rf "$species_dir"
fi
done < "$GENOME_DIR/species_taxids.txt"
echo "Downloaded and organized genomes for $SPECIES"
}
get_salmonella_serotype_taxid() {
esearch -db taxonomy -query "Salmonella enterica subsp. enterica[Subtree]" | efetch -format xml | xtract -pattern Taxon -element TaxId,ScientificName | grep -v -E "str\.|var\." | grep -v "Salmonella enterica subsp. enterica$" > "$GENOME_DIR/salmonella_serotype_taxids.txt"
echo "Saved all serotype Taxon IDs in salmonella_serotype_taxids.txt"
# Remove all Salmonella enterica subsp. enterica serovar prefix
sed -i 's/Salmonella enterica subsp. enterica serovar //g' "$GENOME_DIR/salmonella_serotype_taxids.txt"
# Delete taxids that are already stored in all_monophasic_Typhimurium_taxids.txt
awk 'NR==FNR {taxids[$1]; next} !($1 in taxids)' "$GENOME_DIR/all_monophasic_Typhimurium_taxids.txt" "$GENOME_DIR/salmonella_serotype_taxids.txt" > tmp && mv tmp "$GENOME_DIR/salmonella_serotype_taxids.txt"
} 
download_salmonella_serotype() {
echo "Downloading all Salmonella enterica serotype complete genomes..."
# Create an associative array to store serotypes and their taxids
declare -A serotype_taxids
while IFS=$'\t' read -r taxid serotype; do
clean_serotype=$(echo "$serotype" | sed 's/ /_/g')
serotype_taxids["$clean_serotype"]="$taxid"
done < "$GENOME_DIR/salmonella_serotype_taxids.txt"
# Process all serotypes other than monophasic Typhimurium
for serotype in "${!serotype_taxids[@]}"; do
local serotype_dir="${GENOME_DIR}/Salmonella/Salmonella_enterica/$serotype"
mkdir -p "$serotype_dir"
if [[ -n "$(find "$serotype_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for Salmonella enterica $serotype, skipping."
continue
fi
taxid="${serotype_taxids[$serotype]}"
echo "Downloading genomes for Salmonella $serotype (TaxID: $taxid)..."
ncbi-genome-download bacteria --taxid "$taxid" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$serotype_dir"
if [[ -d "$serotype_dir/genbank" ]]; then
find "$serotype_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$1" && mv "${1%.gz}" "'"$serotype_dir"'"' _ {} \;
rm -rf "$serotype_dir/genbank"
fi
echo "Downloaded genomes for Salmonella $serotype"
done
local monophasic_dir="${GENOME_DIR}/Salmonella/Salmonella_enterica/monophasic_Typhimurium"
mkdir -p "$monophasic_dir"
if [[ -n "$(find "$monophasic_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for Salmonella monophasic Typhimurium, skipping."
else
echo "Downloading all TaxIDs for monophasic Typhimurium"
taxid_list=$(tr '\n' ',' < "$MONOPHASIC_TYPHIMURIUM_TAXIDS" | sed 's/,$//')
ncbi-genome-download bacteria --taxid "$taxid_list" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$serotype_dir"
find "$monophasic_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$1" && mv "${1%.gz}" "'"$monophasic_dir"'"' _ {} \;
rm -rf "$monophasic_dir/genbank"
fi
echo "Downloaded genomes for Salmonella monophasic Typhimurium."
find "${GENOME_DIR}/Salmonella/Salmonella_enterica" -type d -empty -delete
}
# Download genomes 
for GENUS in "${EB_GENUS[@]}"; do
download_genus "$GENUS"
done
get_species_taxids
download_species_by_taxid 
get_salmonella_serotype_taxid
download_salmonella_serotype
