#!/bin/bash
#SBATCH --job-name=build_EB_complete_genomes_db
#SBATCH --output=slurm_build_EB_db_%j.out
#SBATCH --error=slurm_build_EB_db_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
set -euo pipefail
# ============================
# Set working directories and paths
# ============================
WORKDIR=$(pwd)
DATABASE_DIR="$WORKDIR/database"
GENOME_DIR="$DATABASE_DIR/complete_genomes"
BLAST_DB_DIR="$DATABASE_DIR/complete_blast_db"
FAILED_FLAG="WORKDIR/build_EB_db_failed.flag"
MONOPHASIC_TYPHIMURIUM_LIST="$DATABASE_DIR/monophasic_Typhimurium_list.txt"
mkdir -p "$GENOME_DIR"
mkdir -p "$BLAST_DB_DIR"
> "$FAILED_FLAG"
# Load all required functions
source "${WORKDIR}/function.sh"

# ============================
# Dynamic throttling based on available CPU cores
# ============================
if [[ -n "${SLURM_CPUS_ON_NODE:-}" ]]; then
TOTAL_CPUS=$SLURM_CPUS_ON_NODE
else
TOTAL_CPUS=$(nproc)
fi

# Static limits
max_parallel_genus=6
max_parallel_serotypes=6

echo "Detected $TOTAL_CPUS logical cores, parallel limits:"
echo " - Genus: $max_parallel_genus; - Species: $max_parallel_species; -Serotypes: $max_parallel_serotypes"

# ============================
# Step 1: Download genus and species in parallel
# ============================
# List of Enterobacteriaceae Genus used to build the database
GENUS=("Escherichia" "Salmonella" "Shigella" "Klebsiella" "Enterobacter" "Cronobacter" "Citrobacter")
echo "Starting genus-level download and species-parallelization"
for GENUS in "${GENUS[@]}"; do
(
if ! download_genus "$GENUS" "$GENOME_DIR"; then
echo "Failed to download genus: $GENUS" >> "$FAILED_FLAG"
exit 1
fi

if ! get_species_list "$GENUS" "$GENOME_DIR"; then
echo "Failed to get species list for: $GENUS" >> "$FAILED_FLAG"
exit 1
fi

if ! download_species "$GENUS" "$GENOME_DIR"; then
echo "Failed to download species under genus: $GENUS" >> "$FAILED_FLAG"
fi
) &
while (( $(jobs -r | wc -l) >= max_parallel_genus )); do
sleep 1
done
done
# Wait for all background jobs to complete
wait

# ============================
# Step 2: Download Salmonella subspecies (sequential) and serotypes in parallel
# ============================
echo "Starting Salmonella subspecies download"
if ! get_salmonella_subsp_list "$GENOME_DIR"; then
echo "Failed to get Salmonella subspecies list" >> "$FAILED_FLAG"
fi

if ! download_salmonella_subsp "$GENOME_DIR"; then
echo "Failed to download Salmonella subspecies" >> "$FAILED_FLAG"
fi

echo "Starting downloading Salmonella serotypes in parallel"

if ! get_salmonella_serotype_list "$GENOME_DIR"; then
echo "Failed to get Salmonella serotype list" >> "$FAILED_FLAG"
fi

if ! download_salmonella_serotype "$GENOME_DIR"; then
echo "Failed to download Salmonella serotypes" >> "$FAILED_FLAG"
fi

# ============================
# Step 3: Fail-safe check before continuing
# ============================
if [[ -s "$FAILED_FLAG" ]]; then
echo "One or more genome downloads failed. See $FAILED_FLAG for details."
echo "BLAST database will NOT be built."
exit 1
fi

# ============================
# Step 4: Organize unclassified genomes and build BLAST database
# ============================
echo "Organizing unclassified genomes"
move_unclassified_genomes "$GENOME_DIR"

echo "Building BLAST databases"
build_blastdb "$GENOME_DIR" "$BLAST_DB_DIR"

echo "Genome download and BLAST database build completed successfully."
