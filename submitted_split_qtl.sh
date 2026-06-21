#!/bin/bash
set -euo pipefail

# =============================================================================
# submit_split_qtl.sh
# =============================================================================

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

INPUT_DIR="${1%/}"
OUTPUT_DIR="$2"

SPLIT_SCRIPT="/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/split_qtl.sh"
JOB_DIR="/oak/stanford/groups/smontgom/hkrupkin/GxE/Code/slurm_logs"
mkdir -p "${JOB_DIR}"

# --- Derive a unique name from the input directory ---
# e.g., metabolome_chunk_1, proteome, etc.
DIR_LABEL=$(basename "${INPUT_DIR}")

# Unique prefix list per input directory
PREFIX_LIST="${JOB_DIR}/split_qtl_prefixes_${DIR_LABEL}.txt"

# --- Discover unique prefixes ---
ls "${INPUT_DIR}"/*.cis_qtl_pairs.*.csv 2>/dev/null \
    | sed 's/\.[0-9][0-9]*\.csv$//' \
    | sort -u \
    > "${PREFIX_LIST}"

N_JOBS=$(wc -l < "${PREFIX_LIST}")

if [[ ${N_JOBS} -eq 0 ]]; then
    echo "ERROR: No .cis_qtl_pairs.*.csv files found in ${INPUT_DIR}"
    exit 1
fi

echo "Found ${N_JOBS} unique celltype/assay combinations in ${DIR_LABEL}:"
cat "${PREFIX_LIST}" | while read -r p; do basename "$p"; done
echo ""
echo "Prefix list: ${PREFIX_LIST}"
echo "Submitting SLURM array 1-${N_JOBS}..."
echo ""

# --- Submit ---
sbatch << EOF
#!/bin/bash
#SBATCH --job-name=split_${DIR_LABEL}
#SBATCH --time=4:00:00
#SBATCH --account=smontgom
#SBATCH --partition=nih_s10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -e ${JOB_DIR}/split_${DIR_LABEL}_%A_%a.err
#SBATCH -o ${JOB_DIR}/split_${DIR_LABEL}_%A_%a.out
#SBATCH --array=1-${N_JOBS}

set -euo pipefail

PREFIX=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "${PREFIX_LIST}")

echo "Task \${SLURM_ARRAY_TASK_ID}: \${PREFIX}"
bash "${SPLIT_SCRIPT}" "\${PREFIX}" "${OUTPUT_DIR}"
EOF

echo "Submitted: ${DIR_LABEL} (${N_JOBS} tasks)"