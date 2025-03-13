#!/bin/bash
# Initial variables
WORKDIR=$(pwd)
GENOME_DIR="$WORKDIR/database/complete_genomes"
ASSEMBLY_LEVEL="complete"
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
ncbi-genome-download bacteria --genera "$genus" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$genus_dir" --verbose
find "$genus_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$0" && mv "${0%.gz}" "'"$genus_dir"'"' {} \;
rm -rf "$genus_dir/genbank"
echo "Downloaded and organized genomes for $genus"
}
# Build a function to get all species taxid
get_species_taxids() {
for GENUS in "${EB_GENUS[@]}"; do
esearch -db taxonomy -query "$GENUS[Subtree]" | efetch -format xml | xtract -pattern Taxon -element TaxId,ScientificName | awk 'NF == 3 && $2 ~ /^[A-Z][a-z]+$/ && $3 ~ /^[a-z]+$/ {print $1, $2, $3}' >> "$GENOME_DIR/species_taxids.txt"
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
echo "Genomes already exist for $species, skipping donwloading"
continue
fi
echo "Downloading genomes for $species"
ncbi-genome-download bacteria --taxid "$TAXID" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$species_dir"
# Move genomic FASTA files to correct location
find "$species_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$0" && mv "${0%.gz}" "'"$species_dir"'"' {} \;
rm -rf "$species_dir/genbank"
done < "$GENOME_DIR/species_taxids.txt"
find "${GENOME_DIR}/$genus" -type d -empty -delete
echo "Downloaded and organized genomes for $species"
}
get_salmonella_serotype_taxid() {
esearch -db taxonomy -query "Salmonella enterica subsp. enterica[Subtree]" | efetch -format xml | xtract -pattern Taxon -element TaxId,ScientificName | grep -v -E "str\.|var\." | \
grep -v "Salmonella enterica subsp. enterica$" > "$GENOME_DIR/salmonella_serotype_taxids.txt"
echo "Saved all serotype Taxon IDs in salmonella_serotype_taxids.txt"
sed -i 's/Salmonella enterica subsp. enterica serovar //g' "$GENOME_DIR/salmonella_serotype_taxids.txt"
} 
download_salmonella_serotype() {
echo "Downloading all Salmonella enterica serotype complete genomes..."
while read -r taxid serotype; do
local clean_serotype=$(echo "$serotype" | sed 's/ /_/g')
local serotype_dir="${GENOME_DIR}/Salmonella/Salmonella_enterica/${clean_serotype}"
mkdir -p "$serotype_dir"
if [[ -n "$(find "$serotype_dir" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for Salmonella enterica $serotype, skipping."
continue
fi
ncbi-genome-download bacteria --taxid "$taxid" --assembly-level "$ASSEMBLY_LEVEL" --formats fasta --section genbank --output-folder "$serotype_dir"
find "$serotype_dir/genbank" -type f -name "*_genomic.fna.gz" -exec sh -c 'gzip -d "$0" && mv "${0%.gz}" "'"$serotype_dir"'"' {} \;
rm -rf "$serotype_dir/genbank"
find "${GENOME_DIR}/Salmonella/Salmonella_enterica" -type d -empty -delete
echo "Downloaded genomes for Salmonella $serotype"
done < "$GENOME_DIR/salmonella_serotype_taxids.txt"
}
# Download genomes 
for GENUS in "${EB_GENUS[@]}"; do
download_genus "$GENUS"
done
get_species_taxids
download_species_by_taxid 
get_salmonella_serotype_taxid
download_salmonella_serotype
