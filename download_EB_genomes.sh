#!/bin/bash
# Initial variables
WORKDIR=$(pwd)
GENOME_DIR="database/EB_complete_genomes"
ASSEMBLY_LEVEL="complete"
# List of Enterobacteriaceae Genera
EB_GENUS=("Escherichia" "Salmonella" "Shigella" "Klebsiella" "Enterobacter" "Cronobacter" "Citrobacter")
source "$WORKDIR/function.sh"
# Build a function to download genomes using datasets (for genus and species only)
download_genus() {
local genus="$1"
local download_dir="$GENOME_DIR/$genus/aggregated"
if [[ -n "$(find "$download_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes exist for $genus, skipping"
return 0
fi
mkdir -p "$download_dir"
local zip_file="${download_dir}.zip"
echo "Downloading genus: $genus"
datasets download genome taxon "$genus" --assembly-level "$ASSEMBLY_LEVEL" --filename "$zip_file"
unzip -q -o "$zip_file" -d "$download_dir"
rm -f "$zip_file"
find "$download_dir" -type f -name "*_genomic.fna" -exec mv {} "$download_dir" \;
rm -rf "$download_dir/ncbi_dataset"
echo "Downloaded and organized genomes for $genus"
}
download_species() {
local species="$1"
local genus=$(echo "$species" | awk '{print $1}')
local clean_species=$(echo "$species" | sed 's/ /_/g')
local taxon_dir="$GENOME_DIR/$genus/$clean_species"
mkdir -p "$taxon_dir/aggregated"
if [[ -n "$(find "$taxon_dir/aggregated" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for $species, skipping donwloading"
return 0
fi
local zip_file="${taxon_dir}.zip"
echo "Downloading genomes for $species"
datasets download genome taxon "$species" --assembly-level "$ASSEMBLY_LEVEL" --filename "$zip_file"
unzip -q -o "$zip_file" -d "$taxon_dir"
rm -rf "$zip_file"
# Move genomic FASTA files to correct location
find "$taxon_dir" -type f -name "*_genomic.fna" -exec mv {} "$taxon_dir" \;
rm -rf "$taxon_dir/ncbi_dataset"
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
while read TAXID SEROTYPE; do
CLEAN_SEROTYPE=$(clean_taxons "$SEROTYPE")
DOWNLOAD_DIR="${GENOME_DIR}/Salmonella/Salmonella_enterica/${CLEAN_SEROTYPE}"
mkdir -p "${DOWNLOAD_DIR}"
if [[ -n "$(find "$DOWNLOAD_DIR" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for Salmonella enterica $SEROTYPE."
continue
fi
ZIP_FILE="${DOWNLOAD_DIR}.zip"
datasets download genome taxon "$TAXID" --assembly-level "$ASSEMBLY_LEVEL" --filename "$ZIP_FILE"
unzip -q -o "$ZIP_FILE" -d "$DOWNLOAD_DIR"
rm -rf "$ZIP_FILE"
find "$DOWNLOAD_DIR/ncbi_dataset/data/" -type f -name "*_genomic.fna" -exec mv {} "$DOWNLOAD_DIR/" \;
rm -rf "$DOWNLOAD_DIR/ncbi_dataset"
done < salmonella_serotype_taxids.txt
}
# Download genomes 
for GENUS in "${EB_GENUS[@]}"; do
mapfile -t SPECIES_LIST < <(esearch -db taxonomy -query "$GENUS[Subtree]" | efetch -format xml | xtract -pattern Taxon -element ScientificName | grep -E "^$GENUS [a-z]+$")
for SPECIES in "${SPECIES_LIST[@]}"; do
download_genomes "$SPECIES" 
done
done
get_salmonella_serotype_taxid
download_salmonella_serotype
