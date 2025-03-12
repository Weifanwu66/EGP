#!/bin/bash
GENOME_DIR="database/EB_complete_genomes"
BLAST_DB_DIR="database/EB_blast_db"
mkdir -p "${BLAST_DB_DIR}"
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
find "${GENOME_DIR}/${GENUS}" -type f -name "*_genomic.fna" -exec cat {} + > "$GENUS_FASTA"
create_blastdb "$GENUS_FASTA" "$GENUS"
# Species-level databases
SPECIES_LIST=$(find "${GENOME_DIR}/${GENUS}" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;)
for SPECIES in $SPECIES_LIST; do
SPECIES_FASTA="${GENOME_DIR}/${GENUS}/${SPECIES}/${SPECIES}_all_genomes.fna"
find "${GENOME_DIR}/${GENUS}/${SPECIES}" -type f -name "*_genomic.fna" -exec cat {} + > "$SPECIES_FASTA"
create_blastdb "$SPECIES_FASTA" "${SPECIES}"
done
done
echo "Processing serotype level database for Salmonella enterica"
SEROTYPE_DIR="${GENOME_DIR}/Salmonella/Salmonella_enterica"
for DIR in "${SEROTYPE_DIR}"/*; do
SEROTYPE=$(basename "$DIR")
SEROTYPE_FASTA="${DIR}/${SEROTYPE}_genomes.fna"
find "${DIR}" -type f -name "*_genomic.fna" -exec cat {} + > "$SEROTYPE_FASTA"
create_blastdb "$SEROTYPE_FASTA" "Salmonella_${SEROTYPE}"
done

