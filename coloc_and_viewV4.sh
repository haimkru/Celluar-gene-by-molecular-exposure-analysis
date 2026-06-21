#!/bin/bash
#SBATCH --job-name=coloc2V3
#SBATCH --time=24:00:00
#SBATCH --account=smontgom
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH -e coloc2V3%A_%a.err
#SBATCH -o coloc2V3%A_%a.out

set -euo pipefail

# Defaults (override with exported env vars if needed)
R_SCRIPT="/labs/smontgom/grps_smontgom/hkrupkin/GxE/Code/coloc_and_viewV4.R"
ALL_HITS_TSV="${ALL_HITS_TSV:-/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/all_cell_types_sig_hits.tsv}"
TERM_LIST="${TERM_LIST:-/labs/smontgom/grps_smontgom/hkrupkin/GxE/Data/coloc_interaction_with_GWASes/term_list_V2.txt}"
P_INT_THRESH="${P_INT_THRESH:-1e-5}"

# First invocation: build term list and submit array.
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    mkdir -p "$(dirname "$TERM_LIST")"

    awk -F'\t' -v pthresh="$P_INT_THRESH" '
        NR==1 {
            for (i=1; i<=NF; i++) {
                if ($i=="term")        term_col=i
                if ($i=="effect_type") effect_col=i
                if ($i=="pval")        pval_col=i
            }
            if (!term_col || !effect_col || !pval_col) {
                print "Missing required columns: term, effect_type, pval" > "/dev/stderr"
                exit 1
            }
            next
        }
        $effect_col=="interaction_effect" && ($pval_col+0)<=pthresh {
            if (!seen[$term_col]++) print $term_col
        }
    ' "$ALL_HITS_TSV" > "$TERM_LIST"

    n_terms=$(wc -l < "$TERM_LIST")
    if [[ "$n_terms" -lt 1 ]]; then
        echo "No qualifying terms found."
        echo "Checked: $ALL_HITS_TSV"
        exit 1
    fi

    echo "Generated term list: $TERM_LIST"
    echo "Terms found: $n_terms"
    cat "$TERM_LIST"
    echo "Submitting array with $n_terms tasks..."
    sbatch --array=1-"$n_terms" --export=ALL,R_SCRIPT="$R_SCRIPT",TERM_LIST="$TERM_LIST" "$0"
    exit 0
fi

# Array task invocation: run one term.
if [[ ! -s "$TERM_LIST" ]]; then
    echo "TERM_LIST missing or empty: $TERM_LIST"
    exit 1
fi

TERM_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TERM_LIST")
if [[ -z "${TERM_NAME:-}" ]]; then
    echo "No term found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

module load R/4.3.3

echo "starting run | task=${SLURM_ARRAY_TASK_ID} | term=${TERM_NAME}"
Rscript "$R_SCRIPT" "--term=${TERM_NAME}"
echo "done :) | term=${TERM_NAME}"