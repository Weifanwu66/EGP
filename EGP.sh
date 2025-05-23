#!/bin/bash
# Set initial variables and paths
MODE="light"
MIN_COVERAGE=80
MIN_IDENTITY=90
WORK_DIR=$(pwd)
TAXON_FILE=""
GENE_FILE=""
DOWNLOAD_FILE=""
FAILED_FLAG="${WORK_DIR}/build_custom_db_failed.flag"
> "$FAILED_FLAG"
OUTPUT_FILE="${WORK_DIR}/result/gene_summary.csv"
DATABASE_DIR="${WORK_DIR}/database"
BLAST_DB_DIR="${DATABASE_DIR}/complete_blast_db"
CUSTOM_GENOMES_DIR="${DATABASE_DIR}/complete_genomes"
BLAST_RESULT_DIR="${WORK_DIR}/result/complete_blast_results"
FILTERED_BLAST_RESULT_DIR="${WORK_DIR}/result/filtered_complete_blast_results"
DRAFT_GENOMES_DIR="${DATABASE_DIR}/draft_genomes"
DRAFT_BLAST_DB_DIR="${DATABASE_DIR}/draft_blast_db"
DRAFT_BLAST_RESULT_DIR="${WORK_DIR}/result/draft_blast_results"
FILTERED_DRAFT_BLAST_RESULT_DIR="${WORK_DIR}/result/filtered_draft_blast_results"
OVERWRITE=false
# job scheduler variables
runtime=24:00:00; hpcmem=360GB; hpcthreads=72; hpc=F; queue=NA; account=NA
# usage setup
usage() {
    echo "Usage: $0 -g GENE_FILE [-t TAXON_FILE] [-d DOWNLOAD_FILE] [-c COVERAGE] [-i IDENTITY] [-o OUTPUT_FILE] [-p HPC_CLUSTER] [-q QUEUE] [-r RUNTIME] [-m MEMORY] [-C THREADS] [-a ACCOUNT] [-H MODE] [-O OVERWRITE] [-h|--help]"
    echo ""
    echo "Required arguments:"
    echo "-g GENE_FILE      : FASTA file containing target gene sequences."
    echo ""
    echo "Optional arguments:"
    echo "-t TAXON_FILE     : File containing target species (one per line) or a single taxon name."
    echo "-d DOWNLOAD_FILE  : File listing target genera (one per line) to download complete genomes."
    echo "-c COVERAGE       : Minimum genome coverage threshold (default: 80%)."
    echo "-i IDENTITY       : Minimum percentage identity threshold (default: 90%)."
    echo "-o OUTPUT_FILE    : Output filename for results (default: gene_summary.tsv)."
    echo "-p HPC_CLUSTER    : HPC system name (optional; future compatibility)."
    echo "-q QUEUE          : Queue/partition name (e.g., ceres, short, long)."
    echo "-r RUNTIME        : SLURM walltime request (e.g., 04:00:00)."
    echo "-m MEMORY         : Memory request for SLURM job (e.g., 16G)."
    echo "-C THREADS        : Number of CPU cores to request for SLURM."
    echo "-a ACCOUNT        : SLURM account/project (if needed)."
    echo "-H MODE           : Analysis mode ('light' or 'heavy'). Default is light."
    echo "-O OVERWRITE      : Set to true to overwrite previous results (default: false)."
    echo "-h, --help        : Show this help message and exit."
}
# Parse argument (adapted)
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -g) GENE_FILE="$2"; shift 2 ;;
        -t) TAXON_FILE="$2"; shift 2 ;;
        -d) DOWNLOAD_FILE="$2"; shift 2 ;;
        -c) MIN_COVERAGE="$2"; shift 2 ;;
        -i) MIN_IDENTITY="$2"; shift 2 ;;
        -o) OUTPUT_FILE="$2"; shift 2 ;;
        -p) hpc="$2"; shift 2 ;;
        -q) queue="$2"; shift 2 ;;
        -r) runtime="$2"; shift 2 ;;
        -m) hpcmem="$2"; shift 2 ;;
        -C) hpcthreads="$2"; shift 2 ;;
        -a) account="$2"; shift 2 ;;
        -H) MODE="$2"; shift 2 ;;
        -O) OVERWRITE="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Invalid option: $1"; usage; exit 1 ;;
    esac
done

#adapted from GEAbash_v1.0.0; seems to be working as expected
#while getopts ':g:t::c::i::o::p::q::r::m::C::a::H::O::h::' flag; do
#  case "${flag}" in
#    g) GENE_FILE="${OPTARG}" ;;    t) TAXON_FILE="${OPTARG}" ;;    c) MIN_COVERAGE="${OPTARG}" ;;
#    i) MIN_IDENTITY="${OPTARG}" ;;    o) OUTPUT_FILE="${OPTARG}" ;;    p) hpc="${OPTARG}" ;;
#    q) queue="${OPTARG}" ;;    r) runtime="${OPTARG}" ;;    m) hpcmem="${OPTARG}" ;;
#    C) hpcthreads="${OPTARG}" ;;    a) account="${OPTARG}" ;;    H) MODE="${OPTARG}" ;;
#    O) OVERWRITE="${OPTARG}" ;;     
#    -h|--help) usage; exit 0 ;;
#    *) echo "Invalid option: $1"; usage
#       exit 1 ;;    esac; done


while [ $hpc == F ]; do

# detect and utilize slurm or sge job manager. tested on atlas/ceres (slurm) and hank (sge) hpcs and putatively any other systems as requested by users.
# progresses to an expected `find` error line 292
# suspected dependencies: 1) conda (unless some testing/development done independent). e.g. the batch script slurm template in GEAbash calls 
# module load apptainer and I expect this slurm template will ultimately need to call module load miniconda (ceres) or module load miniconda3 (atlas)
# 2) The batch template will need to be in the same directory as EGP.sh and both will need to be in $(pwd)
# more testing to be done after have THE database.
if [ ! $queue == "NA" ]; then 

#detect job submission system
if [ -n "$(sinfo 2>/dev/null)" ]; then submsys=slurm
#detect sge
elif [ -n "$(qhost 2>/dev/null)" ]; then submsys=sge
#exit if queue specified but no job scheduler detected
else echo -n "queue specified but hpc system "
echo "uncertain. exiting"; exit 1; fi

#if sge, remove any logs from previous runs
if [[ " ${submsys} " = " sge " ]]; then 
if [ -f sge2.sh ]; then rm sge2_EGP.* >/dev/null 2>&1; sleep 60; fi
#then go get a new sge template file to use for THIS EGP run
cp sge.sh sge2.sh
#if slurm, remove any logs from previous runs
elif [[ " ${submsys} " = " slurm " ]]; then
if [ -f slurm2.sh ]; then rm slurm2*.out >/dev/null 2>&1; sleep 60; fi
#then go get a new slurm template file to use for THIS GEAbash run
cp slurm.sh slurm2.sh
#otherwise exit. this line should be unnecessary
else echo "system uncertain. exiting"; exit 1
fi

###if sge...
if [[ " ${submsys} " = " sge " ]]
#alert the user
then echo -n "preparing to run EGP in hpc cluster mode. "
echo -n "EGP log outputs will be in the hpc submission system "
echo "log files for sge. e.g., EGP.e* & EGP.o*"
#edit the sge template with local variables and user arguments
sed -i "s/name/sge2_EGP/g" sge2.sh #local
sed -i "s/queue/$queue/g" sge2.sh #user
sed -i "s/runtime/$runtime/g" sge2.sh #user
sed -i "s/RAM/$hpcmem/g" sge2.sh #user
sed -i "s/hpctasks/$hpcthreads/g" sge2.sh #user
#write the EGP command to the sge template
#to use this design syntax, all user options need single dash shortcuts
#e.g --mode heavy and --overwrite become -H heavy and -O true respectively.
#lines 201-205 assume -m <RAM> -q <queue> -r <runtime> in 
#user supplied arguments to EGP or EGP defaults
#lines 214,221 assume -g -c -i -o -H -O -t in 
#user supplied arguments to EGP or EGP defaults
sed -i \
's%command%bash EGP.sh -g "${GENE_FILE}" -c "${MIN_COVERAGE}" -i "${MIN_IDENTITY}" -o "${OUTPUT_FILE}" -H "${MODE}" -O "${OVERWRITE}" -t "${TAXON_FILE}" -d "${DOWNLOAD_FILE}" -p T%g' \
sge2.sh
#write the needed variables to the sge template
#the last variable tells EGP that ${"hpc"} = T
#so that this entire while loop will be skipped
#by the EGP resubmission
sed -i \
"s%vars%OVERWRITE='$OVERWRITE'; GENE_FILE='$GENE_FILE'; MIN_COVERAGE='$MIN_COVERAGE'; MIN_IDENTITY='$MIN_IDENTITY'; OUTPUT_FILE='$OUTPUT_FILE'; MODE='$MODE'; TAXON_FILE='$TAXON_FILE'; DOWNLOAD_FILE='$DOWNLOAD_FILE'%g" \
sge2.sh

#make sure the user provided account is written to the template
if [ $account == "NA" ]
then other='##'
sed -i "s/account/$other/g" sge2.sh
sed -i "s/-P #/###/g" sge2.sh
else sed -i "s/account/$account/g" sge2.sh #an optional user supplied variable
fi
#submit the sge batch job containing the EGP
#command from line 214
qsub sge2.sh
#exit normally without error
exit 0

###if slurm...
else echo -n "preparing to run EGP in hpc cluster mode. "
echo -n "EGP log outputs will be in the hpc submission system "
echo "log files for slurm. e.g., slurm2-*.out"
#edit the slurm template with local variables and user arguments
sed -i "s/name/slurm2_EGP/g" slurm2.sh #local
sed -i "s/queue/$queue/g" slurm2.sh #user
sed -i "s/runtime/$runtime/g" slurm2.sh #user
sed -i "s/RAM/$hpcmem/g" slurm2.sh #user
sed -i "s/hpctasks/$hpcthreads/g" slurm2.sh #user
sed -i 's/other/$other/g' slurm2.sh #local.
#write the GEA command to the slurm template
sed -i \
's%command%bash EGP.sh -g "${GENE_FILE}" -c "${MIN_COVERAGE}" -i "${MIN_IDENTITY}" -o "${OUTPUT_FILE}" -H "${MODE}" -O "${OVERWRITE}" -t "${TAXON_FILE}" -d "${DOWNLOAD_FILE}" -p T%g' \
slurm2.sh
#write the needed variables to the slurm template
#the last variable tells EGP that ${"hpc"} = T
#so that this entire while loop will be skipped
#by the EGP resubmission
sed -i \
"s%vars%OVERWRITE='$OVERWRITE'; GENE_FILE='$GENE_FILE'; MIN_COVERAGE='$MIN_COVERAGE'; MIN_IDENTITY='$MIN_IDENTITY'; OUTPUT_FILE='$OUTPUT_FILE'; MODE='$MODE'; TAXON_FILE='$TAXON_FILE'; DOWNLOAD_FILE='$DOWNLOAD_FILE'%g" \
slurm2.sh

#make sure the user provided account is written to the template
if [ $account == "NA" ]
then other='##'
sed -i "s/account/$other/g" slurm2.sh
sed -i "s/-A #/###/g" slurm2.sh
else sed -i "s/account/$account/g" slurm2.sh
fi
#submit the slurm batch job containing the EGP
#command from line 250
sbatch slurm2.sh
#exit normally without error
exit 0
#finish the if statement started line 195.
fi
#finish the if statement started line 170.
fi
#if the user has NOT specified a queue then prepare
#to run EGP normally
if [ $queue == "NA" ]
then echo "EGP continuing without job submission system"
#if a job submission system is detected,
#alert the user and exit with error.
if [ -n "$(sinfo 2>/dev/null)" ]
then echo -n "queue not specified for running "
echo "in hpc mode. exiting"; exit 1
elif [ -n "$(qhost 2>/dev/null)" ]
then echo -n "queue not specified for running "
echo "in hpc mode. exiting"; exit 1
fi
#finish the if statement started line 161
fi
#if made it this far, set hpc to T to break while loop.
hpc=T
done

if [[ "$OVERWRITE" == true ]]; then
rm -rf "$BLAST_RESULT_DIR" "$FILTERED_BLAST_RESULT_DIR" "$DRAFT_BLAST_RESULT_DIR" "$DRAFT_BLAST_DB_DIR" "$FILTERED_BLAST_RESULT_DIR" "$DRAFT_GENOMES_DIR" "$FILTERED_DRAFT_BLAST_RESULT_DIR" "$OUTPUT_FILE"
fi
# Ensure gene file is provided
if [[ -z "$GENE_FILE" ]]; then
echo "Error! Please provide a gene sequence file!"
usage
exit 1
fi
# Ensure if heavy mode is enabled, a taxon file must be provided
if [[ "$MODE" == "heavy" && -z "$TAXON_FILE" ]]; then
echo "Error: A taxon file is required in heavy mode."
exit 1
fi
# Output directory setup
mkdir -p "$BLAST_RESULT_DIR" "$FILTERED_BLAST_RESULT_DIR"
source "${WORK_DIR}/function.sh" || { echo "Error sourcing function.sh". exit 1; }
# Handle custom download if provided
if [[ -n "$DOWNLOAD_FILE" ]]; then
echo "Custom panel download requested."
GENOME_DIR="$CUSTOM_GENOMES_DIR"
mkdir -p "$GENOME_DIR"
mkdir -p "$BLAST_DB_DIR"
max_parallel_jobs=8
while IFS= read -r raw_line || [ -n "$raw_line" ]; do
taxon=$(echo "$raw_line" | tr -d '\r')
[[ -z "$taxon" ]] && continue
(
if [[ "$taxon" == "Salmonella" || "$taxon" == "Salmonella enterica" ]]; then
echo "$taxon detected. Downloading all subspecies and serotypes."
if ! get_salmonella_subsp_list "$GENOME_DIR"; then
echo "Failed to get Salmonella subspecies list" >> "$FAILED_FLAG"; exit 1
fi
if ! download_salmonella_subsp "$GENOME_DIR"; then
echo "Failed to download Salmonella subspecies" >> "$FAILED_FLAG"; exit 1
fi
if ! get_salmonella_serotype_list "$GENOME_DIR"; then
echo "Failed to get Salmonella serotype list" >> "$FAILED_FLAG"; exit 1
fi
if ! download_salmonella_serotype "$GENOME_DIR"; then
echo "Failed to download Salmonella serotypes" >> "$FAILED_FLAG"; exit 1
fi
elif [[ "$taxon" == "Salmonella enterica subsp. enterica" ]]; then
echo "Downloading all serotypes under 'Salmonella enterica subsp. enterica'"
if ! get_salmonella_serotype_list "$GENOME_DIR"; then
echo "Failed to get serotype list" >> "$FAILED_FLAG"; exit 1
fi
if ! download_salmonella_serotype "$GENOME_DIR"; then
echo "Failed to download serotypes" >> "$FAILED_FLAG"; exit 1
fi
elif [[ "$taxon" =~ ^Salmonella( enterica subsp\.?\ enterica serovar)? ([A-Z][a-zA-Z0-9_]+)$ ]]; then
serotype="${BASH_REMATCH[2]}"
echo "Downloading $serotype"
echo "$serotype" > "$GENOME_DIR/temp_serotype_list.txt"
if ! download_salmonella_serotype "$GENOME_DIR" "$GENOME_DIR/temp_serotype_list.txt"; then
echo "Failed to download serotype: $serotype" >> "$FAILED_FLAG"
fi
rm -f "$GENOME_DIR/temp_serotype_list.txt"
elif [[ "$taxon" =~ ^[A-Z][a-z]+$ ]]; then
GENUS_DIR="${GENOME_DIR}/${taxon}"
if [[ -d "$GENUS_DIR" ]]; then
echo "Genus $taxon already downloaded. Skipping."
else
if ! download_genus "$taxon" "$GENOME_DIR"; then
echo "Failed to download genus: $taxon" >> "$FAILED_FLAG"; exit 1
fi
if ! get_species_list "$taxon" "$GENOME_DIR"; then
echo "Failed to get species list: $taxon" >> "$FAILED_FLAG"; exit 1
fi
if ! download_species "$GENOME_DIR/species_list.txt" "$GENOME_DIR"; then
echo "Failed to download species for $taxon" >> "$FAILED_FLAG"
fi
fi
elif [[ "$taxon" =~ ^[A-Z][a-z]+\ [a-z]+$ ]]; then
SPECIES_DIR="${GENOME_DIR}/${taxon// /_}"
if [[ -d "$SPECIES_DIR" ]]; then
echo "Species $taxon already downloaded. Skipping."
else
if ! download_species "$taxon" "$GENOME_DIR"; then
echo "Failed to download species: $taxon" >> "$FAILED_FLAG"
fi
fi
fi
) &
while (( $(jobs -r | wc -l) >= max_parallel_jobs )); do
sleep 1
done
done < "$DOWNLOAD_FILE"
wait
echo "Organizing unclassified genomes"
move_unclassified_genomes "$GENOME_DIR"
echo "Building BLAST databases for custom panel"
build_blastdb "$GENOME_DIR" "$BLAST_DB_DIR"
if [[ -s "$FAILED_FLAG" ]]; then
echo "Custom panel completed with some failures. See $FAILED_FLAG"
else
echo "Custom panel build complete."
fi
else
GENOME_DIR=""
echo "Default mode: using prebuilt BLAST database."
fi
# Set delimiter as a space
if [[ -n "$TAXON_FILE" ]]; then
if [[ ! -f "$TAXON_FILE" ]]; then
tmp_taxon_file=$(mktemp)
echo "$TAXON_FILE" > "$tmp_taxon_file"
TAXON_FILE="$tmp_taxon_file"
trap "rm -f '$tmp_taxon_file'" EXIT
fi
DELIMITER=" "
sed 's/[\t,]\+/ /g' "$TAXON_FILE" > "${TAXON_FILE}_processed"
# Ensure if a taxon file is provided, the genus should be one of the target Enterobacteriaceae
if [[ -z "$DOWNLOAD_FILE" ]]; then
awk -v delim="$DELIMITER" '
BEGIN { FS=delim; OFS=delim }
{ 
gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
gsub(/^[[:space:]]+|[[:space:]]+$/, "", $2)
genus = $1
species_or_serotype = ($2 != "") ? $2 : ""
if (genus !~ /^(Salmonella|Escherichia|Citrobacter|Enterobacter|Klebsiella|Shigella|Cronobacter)$/) {
print "Error: Invalid genus in taxon file. Allowed: Salmonella, Escherichia, Citrobacter, Enterobacter, Klebsiella, Shigella, Cronobacter in default. To proceed with other genus, use custom panel flag -d."
print "Your line:", $0
exit 1 
}
}' "${TAXON_FILE}_processed" || exit 1
else
echo "Custom database mode detected - skipping taxon genus restriction check."
fi
rm -f "${TAXON_FILE}_processed"
fi
# start with core processing functions
process_complete_genomes() {
local taxon="$1"
echo "Processing complete genomes for $taxon"
local clean_taxon="$(extract_taxon_info "$taxon")"
local blast_output_name="${clean_taxon// /_}"
blast_output_name="${blast_output_name//./}"
local blast_output="${BLAST_RESULT_DIR}/${blast_output_name}_complete_blast_results.txt"
if [[ ! -s "$blast_output" ]]; then
perform_blast "$GENE_FILE" "$MIN_IDENTITY" "$BLAST_RESULT_DIR" "" "$taxon"
echo "BLAST result saved to $blast_output"
else
echo "Skipping BLAST: results already exist for $taxon"
fi
for blast_result_file in "$blast_output"; do
filter_blast_results "$blast_result_file" "$FILTERED_BLAST_RESULT_DIR" "$MIN_COVERAGE" "complete"
done
echo "Finished processing complete genomes for $taxon"
}

process_draft_genomes() {
local taxon="$1"
local clean_taxon="$(extract_taxon_info "$taxon")"
local blast_db_name="${clean_taxon// /_}"
blast_db_name="${blast_db_name//./}"
local total_draft_genomes=$(get_total_genomes_count "$taxon" "contig")
read -r sample_size iterations <<< "$(calculate_sample_size_and_iterations "$total_draft_genomes")"
echo "Processing $taxon | Total draft genomes: $total_draft_genomes. Running $iterations iterations (max 20)."
for ((i=1; i<=iterations; i++)); do
echo "Starting iterations $i/$iterations for $taxon"
download_random_draft_genomes "$taxon" "$sample_size" "$DRAFT_GENOMES_DIR" "$i"
perform_blast "$GENE_FILE" "$MIN_IDENTITY" "$DRAFT_BLAST_RESULT_DIR" "$i" "$taxon"
mkdir -p "${DRAFT_BLAST_RESULT_DIR}"
local blast_result_file="${DRAFT_BLAST_RESULT_DIR}/${blast_db_name}/iteration_${i}_draft_blast_results.txt"
filter_blast_results "$blast_result_file" "$FILTERED_DRAFT_BLAST_RESULT_DIR" "$MIN_COVERAGE" "draft"
done
}

# Set up the output file headers
if [[ "$MODE" == "heavy" ]]; then
echo -e "Organism,Gene_ID,Min_percentage_of_coverage,Min_percentage_of_identity,Total_draft_genomes,Total_complete_genomes,Draft_genomes_sample_size,Number_of_iterations,Complete_genomes_with_target_genes,Draft_genomes_with_target_genes,Percentage_with_target_genes_complete_genomes,Percentage_with_target_genes_draft_genomes" > "$OUTPUT_FILE"
else
echo -e "Organism,Gene_ID,Min_percentage_of_coverage,Min_percentage_of_identity,Total_draft_genomes,Total_complete_genomes,Complete_genomes_with_target_genes,Percentage_with_target_genes_complete_genomes" > "$OUTPUT_FILE"
fi
# Initialize TAXON_LIST
if [[ -n "$TAXON_FILE" ]]; then
TAXON_LIST=$(< "$TAXON_FILE")
else
# Rmove numeric suffix and convert underscores and dots to spaces
TAXON_LIST=$(find "$BLAST_DB_DIR" -name "*.nsq" -exec basename {} .nsq \; | sed -E '/[._][0-9]{2}$/s/[._][0-9]{2}$//' | sort -u |  sed 's/_/ /g' | sed -E 's/subsp /subsp. /g')
fi
# Process each taxon in TAXON_LIST
while IFS= read -r taxon; do
[[ -z "$taxon" ]] && continue
echo "Processing $taxon in $MODE mode"
process_complete_genomes "$taxon"
[[ "$MODE" == "heavy" ]] && process_draft_genomes "$taxon"
TOTAL_COMPLETE_GENOMES=$(get_total_genomes_count "$taxon" "complete") 
TOTAL_DRAFT_GENOMES=$(get_total_genomes_count "$taxon" "contig")
clean_taxon="$(extract_taxon_info "$taxon")"
local_taxon="${clean_taxon// /_}"
local_taxon="${local_taxon//./}"
if [[ "$MODE" == "heavy" && "$TOTAL_DRAFT_GENOMES" -gt 0 ]]; then
read -r DRAFT_SAMPLE_SIZE ITERATIONS <<< "$(calculate_sample_size_and_iterations "$TOTAL_DRAFT_GENOMES")"
else
DRAFT_SAMPLE_SIZE=0
ITERATIONS=0
fi
GENE_WITH_HITS=$(awk '{print $1}' "$FILTERED_BLAST_RESULT_DIR/filtered_${local_taxon}_complete_blast_results.txt" 2>/dev/null)
[[ "$MODE" == "heavy" ]] && GENE_WITH_HITS+=$'\n'$(awk '{print $1}' "$FILTERED_DRAFT_BLAST_RESULT_DIR/$local_taxon"/* 2>/dev/null)
GENE_WITH_HITS=$(echo "$GENE_WITH_HITS" | sort -u)
# Process all genes in query gene file
mapfile -t ALL_GENES < <(grep "^>" "$GENE_FILE" | sed 's/>//' | awk '{print $1}')
for GENE_ID in "${ALL_GENES[@]}"; do
if grep -qx "$GENE_ID" <<< "$GENE_WITH_HITS"; then
COMPLETE_GENOMES_WITH_TARGET_GENES=$(awk -v gene="$GENE_ID" '$1 == gene {print $2}' "$FILTERED_BLAST_RESULT_DIR/filtered_${local_taxon}_complete_blast_results.txt" 2>/dev/null | sort -u | wc -l)
if [[ "$MODE" == "heavy" && "$TOTAL_DRAFT_GENOMES" -gt 0 ]]; then
TOTAL_DRAFT_GENOMES_WITH_TARGET_GENES=0
for ((i=1; i<="$ITERATIONS"; i++)); do
DRAFT_GENOMES_WITH_TARGET_GENES=$(awk -v gene="$GENE_ID" '$1 == gene {print $2}' "$FILTERED_DRAFT_BLAST_RESULT_DIR/${local_taxon}/filtered_iteration_${i}_draft_blast_results.txt" 2>/dev/null | cut -c1-8 | sort -u | wc -l)
TOTAL_DRAFT_GENOMES_WITH_TARGET_GENES=$((TOTAL_DRAFT_GENOMES_WITH_TARGET_GENES + DRAFT_GENOMES_WITH_TARGET_GENES))
done
AVERAGE_DRAFT_GENOMES_WITH_TARGET_GENES=$(echo "$TOTAL_DRAFT_GENOMES_WITH_TARGET_GENES/$ITERATIONS" | bc 2>/dev/null)
PERCENT_WITH_TARGET_GENES_DRAFT_GENOMES=$(echo "scale=2; ($AVERAGE_DRAFT_GENOMES_WITH_TARGET_GENES/$DRAFT_SAMPLE_SIZE) * 100" | bc 2>/dev/null)
else
DRAFT_GENOMES_WITH_TARGET_GENES=0
PERCENT_WITH_TARGET_GENES_DRAFT_GENOMES=0
fi
if [[ "$TOTAL_COMPLETE_GENOMES" -gt 0 ]]; then
PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES=$(echo "scale=2; ($COMPLETE_GENOMES_WITH_TARGET_GENES/$TOTAL_COMPLETE_GENOMES) * 100" | bc 2>/dev/null)
else
PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES=0
fi
if [[ "$MODE" == "heavy" ]]; then
echo -e "$taxon,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$DRAFT_SAMPLE_SIZE,$ITERATIONS,$COMPLETE_GENOMES_WITH_TARGET_GENES,$AVERAGE_DRAFT_GENOMES_WITH_TARGET_GENES,${PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES}%,${PERCENT_WITH_TARGET_GENES_DRAFT_GENOMES}%" >> "$OUTPUT_FILE"
else
echo -e "$taxon,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$COMPLETE_GENOMES_WITH_TARGET_GENES,${PERCENT_WITH_TARGET_GENES_COMPLETE_GENOMES}%"  >> "$OUTPUT_FILE"
fi
else
# If no hits were found for the gene in a serotype, output is recorded as 0%
if [[ "$MODE" == "heavy" ]]; then
echo -e "$taxon,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,$DRAFT_SAMPLE_SIZE,$ITERATIONS,0,0,0%,0%" >> "$OUTPUT_FILE"
else
echo -e "$taxon,$GENE_ID,$MIN_COVERAGE,$MIN_IDENTITY,$TOTAL_DRAFT_GENOMES,$TOTAL_COMPLETE_GENOMES,0,0%" >> "$OUTPUT_FILE"
fi
fi
done
done <<< "$TAXON_LIST"
echo "Analysis complete. Results saved in $OUTPUT_FILE"
