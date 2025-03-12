#!/bin/bash
# Initial variables
WORKDIR=$(pwd)
GENOME_DIR="database/EB_complete_genomes"
ASSEMBLY_LEVEL="complete"
# List of Enterobacteriaceae Genera
EB_GENUS=("Escherichia" "Salmonella" "Shigella" "Klebsiella" "Enterobacter" "Cronobacter" "Citrobacter")
source "$WORKDIR/function.sh"
# Build a function to download genus-level genomes using datasets
download_genus() {
local genus="$1"
local genus_dir="$GENOME_DIR/$genus/aggregated"
if [[ -n "$(find "$genus_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes exist for $genus, skipping"
return
fi
mkdir -p "$genus_dir"
local zip_file="${genus_dir}/${genus}.zip"
echo "Downloading genus: $genus"
datasets download genome taxon "$genus" --assembly-level "$ASSEMBLY_LEVEL" --filename "$zip_file"
unzip -q -o "$zip_file" -d "$genus_dir"
rm -f "$zip_file"
find "$genus_dir/ncbi_dataset/data/" -type f -name "*_genomic.fna" -exec mv {} "$genus_dir" \;
rm -rf "$genus_dir/ncbi_dataset"
echo "Downloaded and organized genomes for $genus"
}
# Build a function to download species-level genomes using datasets
download_species() {
local species="$1"
local genus=$(echo "$species" | awk '{print $1}')
local clean_species=$(echo "$species" | sed 's/ /_/g')
local species_dir="$GENOME_DIR/$genus/$clean_species"
[[ "$species" == "Salmonella enterica" ]] && species_dir="$species_dir/aggregated"
mkdir -p "$species_dir"
if [[ -n "$(find "$species_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for $species, skipping donwloading"
return
fi
local zip_file="${species_dir}/${clean_species}.zip"
echo "Downloading genomes for $species"
datasets download genome taxon "$species" --assembly-level "$ASSEMBLY_LEVEL" --filename "$zip_file"
unzip -q -o "$zip_file" -d "$species_dir"
rm -rf "$zip_file"
# Move genomic FASTA files to correct location
find "$species_dir/ncbi_dataset/data/" -type f -name "*_genomic.fna" -exec mv {} "$species_dir" \;
rm -rf "$species_dir/ncbi_dataset"
find "$species_dir" -type d -empty -delete
echo "Downloaded and organized genomes for $species"
}
get_salmonella_serotype_taxid() {
esearch -db taxonomy -query "Salmonella enterica subsp. enterica[Subtree]" | efetch -format xml | xtract -pattern Taxon -element TaxId,ScientificName | grep -v -E "str\.|var\." | \
grep -v "Salmonella enterica subsp. enterica$" > salmonella_serotype_taxids.txt
echo "Saved all serotype Taxon IDs in salmonella_serotype_taxids.txt"
sed -i 's/Salmonella enterica subsp. enterica serovar //g' salmonella_serotype_taxids.txt
} 
download_salmonella_serotype() {
echo "Downloading all Salmonella enterica serotype complete genomes..."
while read -r taxid serotype; do
local clean_serotype=$(clean_taxons "$serotype")
local serotype_dir="${GENOME_DIR}/Salmonella/Salmonella_enterica/${CLEAN_SEROTYPE}"
mkdir -p "$serotype_dir"
if [[ -n "$(find "$serotype_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for Salmonella enterica $serotype, skipping."
continue
fi
zip_file="${serotype_dir}.zip"
datasets download genome taxon "$taxid" --assembly-level "$ASSEMBLY_LEVEL" --filename "$zip_file"
unzip -q -o "$zip_file" -d "$serotype_dir"
rm -rf "$zip_file"
find "$serotype_dir/ncbi_dataset/data/" -type f -name "*_genomic.fna" -exec mv {} "$serotype_dir" \;
rm -rf "$serotype_dir/ncbi_dataset"
find "$serotype_dir" -type d -empty -delete
echo "Downloaded genomes for Salmonella $serotype"
done < salmonella_serotype_taxids.txt
}
# Download genomes 
for GENUS in "${EB_GENUS[@]}"; do
download_genus "$GENUS"
mapfile -t SPECIES_LIST < <(esearch -db taxonomy -query "$GENUS[Subtree]" | efetch -format xml | xtract -pattern Taxon -element ScientificName | grep -E "^$GENUS [a-z]+$")
for SPECIES in "${SPECIES_LIST[@]}"; do
download_species "$SPECIES" 
done
done
get_salmonella_serotype_taxid
download_salmonella_serotype
