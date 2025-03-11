#!/bin/bash
# Initial variables
WORKDIR=$(pwd)
GENOME_DIR="database/EB_complete_genomes"
ASSEMBLY_LEVEL="complete"
# List of Enterobacteriaceae Genera
EB_GENUS=("Escherichia" "Salmonella" "Shigella" "Klebsiella" "Enterobacter" "Cronobacter" "Citrobacter")
source "$WORKDIR/function.sh"
# Build a function to download genomes using datasets (for genus and species only)
download_genomes() {
TAXON="$1"
CLEAN_TAXON=$(clean_taxons "$TAXON" )
#Define download directory based on level
GENUS=$(echo "$TAXON" | awk '{print $1}')
SPECIES=$(echo "$TAXON" | awk '{print $1, $2}')
CLEAN_SPECIES=$(clean_taxons "$SPECIES")
if [[ "$GENUS" == "Salmonella" && "$SPECIES" == "Salmonella enterica" ]]; then
return
elif [[ "$GENUS" == "Salmonella" && "$SPECIES" == "Salmonella bongori" ]]; then
DOWNLOAD_DIR="${GENOME_DIR}/Salmonella/Salmonella_bongori"
else
DOWNLOAD_DIR="${GENOME_DIR}/${GENUS}/${CLEAN_SPECIES}"
fi
if [[ -n "$(find "$DOWNLOAD_DIR" -maxdepth 1 -type f -name "*_genomic.fna" 2>/dev/null)" ]]; then
echo "Genomes already exist for $TAXON, skipping donwloading"
return
fi
ZIP_FILE="${DOWNLOAD_DIR}.zip"
mkdir -p "${DOWNLOAD_DIR}"
echo "Downloading genomes for $TAXON..."
datasets download genome taxon "$TAXON" --assembly-level "$ASSEMBLY_LEVEL" --filename "${ZIP_FILE}"
unzip -q -o "${ZIP_FILE}" -d "${DOWNLOAD_DIR}"
rm -rf "${ZIP_FILE}"
# Move genomic FASTA files to correct location
find "${DOWNLOAD_DIR}" -type f -name "*_genomic.fna" -exec mv {} "${DOWNLOAD_DIR}" \;
rm -rf "${DOWNLOAD_DIR}/ncbi_dataset"
echo "Downloaded and organized genomes for $TAXON"
}
get_salmonella_serotype_taxid() {
esearch -db taxonomy -query "Salmonella enterica subsp. enterica[Subtree]" | efetch -format xml | xtract -pattern Taxon -element TaxId,ScientificName | grep -v -E "str\.|var\."  > salmonella_serotype_taxids.txt
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
#download_genomes "Salmonella bongori"
#get_salmonella_serotype_taxid
download_salmonella_serotype
