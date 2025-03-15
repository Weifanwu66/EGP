#!/bin/bash
WORK_DIR=$(pwd)
GENOME_DIR="$WORK_DIR/database/complete_genomes"
BLAST_DB_DIR="$WORK_DIR/database/complete_blast_db"
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
GENUS_FASTA="${GENOME_DIR}/${GENUS}/aggregated/${GENUS}_all_genomes.fna"
find "${GENOME_DIR}/${GENUS}/aggregated" -type f -name "*_genomic.fna" -exec cat {} + > "$GENUS_FASTA"
create_blastdb "$GENUS_FASTA" "$GENUS"
# Species-level databases
SPECIES_LIST=$(find "${GENOME_DIR}/${GENUS}" -mindepth 1 -maxdepth 1 -type d \( ! -name "aggregated" \) -exec basename {} \;)
for SPECIES in $SPECIES_LIST; do
SPECIES_DIR="${GENOME_DIR}/${GENUS}/${SPECIES}"
if [[ "$SPECIES" == "Salmonella_enterica" ]]; then
AGGREGATED_FASTA="${SPECIES_DIR}/aggregated/Salmonella_enterica_all_genomes.fna"
find "${SPECIES_DIR}/aggregated" -type f -name "*_genomic.fna" -exec cat {} + > "$AGGREGATED_FASTA"
create_blastdb "$AGGREGATED_FASTA" "Salmonella_enterica"
for SEROTYPE_DIR in "${SPECIES_DIR}"/*/; do
SEROTYPE=$(basename "$SEROTYPE_DIR")
if [[ "$SEROTYPE" != "aggregated" ]]; then
SEROTYPE_FASTA="$SEROTYPE_DIR/${SEROTYPE}_genomes.fna"
find "$SEROTYPE_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$SEROTYPE_FASTA"
create_blastdb "$SEROTYPE_FASTA" "Salmonella_${SEROTYPE}"
fi
done
else
SPECIES_FASTA="${SPECIES_DIR}/${SPECIES}_all_genomes.fna"
find "$SPECIES_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$SPECIES_FASTA"
create_blastdb "$SPECIES_FASTA" "$SPECIES"
fi
done
done

