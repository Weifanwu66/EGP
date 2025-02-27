#!/bin/bash
GENOME_DIR="ncbi_genomes"
BLAST_DB_DIR="ncbi_blast_db"
EB_GENUS=("Escherichia" "Salmonella" "Shigella" "Klebsiella" "Enterobacter" "Cronobacter" "Citrobacter")
create_blastdb() {
local INPUT_FASTA="$1"
local DB_NAME="$2"
echo "Creating BLAST database for $DB_NAME..."
makeblastdb -in "$INPUT_FASTA" -dbtype nucl -out "${BLAST_DB_DIR}/${DB_NAME}"
echo "BLAST database created for $DB_NAME at ${BLAST_DB_DIR}/${DB_NAME}"
}
# Create BLAST database for all genus
for GENUS in "${EB_GENUS[@]}"; do
echo "Processing genus: $GENUS"
GENUS_FASTA="${GENOME_DIR}/${GENUS}/${GENUS}_all_genomes.fna"
cat "${GENOME_DIR}/${GENUS}"/*/*_genomic.fna > "$GENUS_FASTA"
create_blastdb "$GENUS_FASTA" "$GENUS"
# Species-level databases
SPECIES_LIST=$(find "${GENOME_DIR}/${GENUS}" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;)
for SPECIES in $SPECIES_LIST; do
SPECIES_FASTA="${GENOME_DIR}/${GENUS}/${SPECIES}/${SPECIES}_all_genomes.fna"
cat "${GENOME_DIR}/${GENUS}/${SPECIES}"/*_genomic.fna > "$SPECIES_FASTA"
create_blastdb "$SPECIES_FASTA" "${SPECIES}"
done
done
echo "Processing serotype level database for Salmonella"
SEROTYPE_DIR="${GENOME_DIR}/Salmonella/Salmonella_enterica"
for DIR in "${SEROTYPE_DIR}"/*; do
SEROTYPE=$(basename "$DIR")
SEROTYPE_FASTA="${SEROTYPE_DIR}/${SEROTYPE}_genomes.fna"
cat "${DIR}"/*_genomic.fna > "$SEROTYPE_FASTA"
create_blastdb "$SEROTYPE_FASTA" "Salmonella_enterica_${SEROTYPE}"
done

