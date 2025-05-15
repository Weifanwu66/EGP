#!/bin/bash
WORKDIR=$(pwd)
DATABASE_DIR="$WORKDIR/database"
GENOME_DIR="$DATABASE_DIR/complete_genomes"
BLAST_DB_DIR="$DATABASE_DIR/complete_blast_db"
MONOPHASIC_TYPHIMURIUM_LIST="$DATABASE_DIR/monophasic_Typhimurium_list.txt"
mkdir -p "$GENOME_DIR"
mkdir -p "$BLAST_DB_DIR"
source "${WORKDIR}/function.sh"
# List of Enterobacteriaceae Genus used to build the database
GENUS=("Escherichia" "Salmonella" "Shigella" "Klebsiella" "Enterobacter" "Cronobacter" "Citrobacter")
echo "Downloading and organizing complete genomes"
max_parallel_jobs=8
# --- Download genomes in parallel ---
for GENUS in "${GENUS[@]}"; do
(
download_genus "$GENUS" "$GENOME_DIR"
get_species_list "$GENUS" "$GENOME_DIR"
download_species "$GENUS" "$GENOME_DIR"
) &
# Limit parallel jobs
while (( $(jobs -r | wc -l) >= max_parallel_jobs )); do
sleep 1
done
done
# Wait for all background jobs to complete
wait

# --- download Salmonella subspecies and serotypes (sequential) ---
get_salmonella_subsp_list "$GENOME_DIR"
download_salmonella_subsp "$GENOME_DIR"
get_salmonella_serotype_list "$GENOME_DIR"
download_salmonella_serotype "$GENOME_DIR"

# --- Organize unclassified genomes ---
echo "Organizing unclassified genomes"
move_unclassified_genomes "$GENOME_DIR"

# --- Build BLAST Databases ---
echo "Building BLAST databases"
build_blastdb "$GENOME_DIR" "$BLAST_DB_DIR"
