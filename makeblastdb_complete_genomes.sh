#!/bin/bash
WORK_DIR=$(pwd)
GENOME_DIR="$WORK_DIR/database/complete_genomes"
BLAST_DB_DIR="$WORK_DIR/database/complete_blast_db"
mkdir -p "${BLAST_DB_DIR}"
create_blastdb() {
local INPUT_FASTA="$1"
local DB_NAME="$2"
echo "Creating BLAST database for $DB_NAME..."
makeblastdb -in "$INPUT_FASTA" -dbtype nucl -out "${BLAST_DB_DIR}/${DB_NAME}"
echo "BLAST database created for $DB_NAME at ${BLAST_DB_DIR}/${DB_NAME}"
}
# Create BLAST database for all genus
find "$GENOME_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r GENUS_DIR; do
GENUS=$(basename "$GENUS_DIR")
echo "Processing genus: $GENUS"
GENUS_FASTA="${GENUS_DIR}/${GENUS}_all_genomes.fna"
find "$GENUS_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$GENUS_FASTA"
#create_blastdb "$GENUS_FASTA" "$GENUS"
# Species-level databases
find "$GENUS_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r SPECIES_DIR; do
SPECIES=$(basename "$SPECIES_DIR")
[[ "$SPECIES" == "unclassified" ]] && continue
SPECIES_FASTA="${SPECIES_DIR}/${SPECIES}_all_genomes.fna"
find "$SPECIES_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$SPECIES_FASTA"
create_blastdb "$SPECIES_FASTA" "$SPECIES"
if [[ "$SPECIES" == "Salmonella_enterica" ]]; then
find "$SPECIES_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r SUBSPECIES_DIR; do
SUBSPECIES=$(basename "$SUBSPECIES_DIR")
[[ "$SUBSPECIES" == "unclassified" ]] && continue
SUBSPECIES_FASTA="${SUBSPECIES_DIR}/Salmonella_enterica_subsp_${SUBSPECIES}_all_genomes.fna"
find "$SUBSPECIES_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$SUBSPECIES_FASTA"
#create_blastdb "$SUBSPECIES_FASTA" "Salmonella_enterica_subsp_${SUBSPECIES}"
if [[ "$SUBSPECIES" == "enterica" ]]; then
find "$SUBSPECIES_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r SEROTYPE_DIR; do
SEROTYPE=$(basename "$SEROTYPE_DIR")
[[ "$SEROTYPE" == "unclassified" ]] && continue
SEROTYPE_FASTA="$SEROTYPE_DIR/${SEROTYPE}_genomes.fna"
find "$SEROTYPE_DIR" -type f -name "*_genomic.fna" -exec cat {} + > "$SEROTYPE_FASTA"
create_blastdb "$SEROTYPE_FASTA" "Salmonella_${SEROTYPE}"
done
fi
done
fi
done
done
