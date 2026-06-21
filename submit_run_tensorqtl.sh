#!/usr/bin/env bash

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
set -o xtrace -o nounset -o errexit

input_directory="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/standardized_phenotypesV2/"
phenotype_list="$(dirname ${input_directory})/phenotype_list"
genotype="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/genotypes/mesa_1331samples.maf01.biallelic.intersect"

find -L "${input_directory}" -type f -name "*.bed" | \
    sort -u > ${phenotype_list}
array_length=$(wc -l ${phenotype_list} | cut -d ' ' -f 1)

code_directory=$(realpath $(dirname ${BASH_SOURCE[0]}))

output_directory="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_outputs/metabolome_chunk_1"
mkdir -p "${output_directory}/logs"
sbatch --output "${output_directory}/logs/%A_%a.log" \
    --error "${output_directory}/logs/%A_%a.log" \
    --array "1-${array_length}" \
    --time 12:00:00 \
    --account sjaiswal \
    --partition gpu \
    --cpus-per-task 1 \
    --gres gpu:1 \
    --mem 64G \
    --job-name tensorqtl \
    "${code_directory}/run_tensorqtl.sh" \
        --genotype "${genotype}" \
        --phenotype_list ${phenotype_list} \
        --interactions "/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/interaction_termsV2/Eigen_Metabolite_Expression_interaction_chunk_1.txt" \
        --output_directory ${output_directory}

output_directory="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_outputs/metabolome_chunk_2"
mkdir -p "${output_directory}/logs"
sbatch --output "${output_directory}/logs/%A_%a.log" \
    --error "${output_directory}/logs/%A_%a.log" \
    --array "1-${array_length}" \
    --time 12:00:00 \
    --account sjaiswal \
    --partition gpu \
    --cpus-per-task 1 \
    --gres gpu:1 \
    --mem 64G \
    --job-name tensorqtl \
    "${code_directory}/run_tensorqtl.sh" \
        --genotype "${genotype}" \
        --phenotype_list ${phenotype_list} \
        --interactions "/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/interaction_termsV2/Eigen_Metabolite_Expression_interaction_chunk_2.txt" \
        --output_directory ${output_directory}

output_directory="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_outputs/metabolome_chunk_3"
mkdir -p "${output_directory}/logs"
sbatch --output "${output_directory}/logs/%A_%a.log" \
    --error "${output_directory}/logs/%A_%a.log" \
    --array "1-${array_length}" \
    --time 12:00:00 \
    --account sjaiswal \
    --partition gpu \
    --cpus-per-task 1 \
    --gres gpu:1 \
    --mem 64G \
    --job-name tensorqtl \
    "${code_directory}/run_tensorqtl.sh" \
        --genotype "${genotype}" \
        --phenotype_list ${phenotype_list} \
        --interactions "/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/interaction_termsV2/Eigen_Metabolite_Expression_interaction_chunk3.txt" \
        --output_directory ${output_directory}

# output_directory="/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_outputs/proteome"
# mkdir -p "${output_directory}/logs"
# sbatch --output "${output_directory}/logs/%A_%a.log" \
#     --error "${output_directory}/logs/%A_%a.log" \
#     --array "1-${array_length}" \
#     --time 12:00:00 \
#     --account sjaiswal \
#     --partition gpu \
#     --cpus-per-task 1 \
#     --gres gpu:1 \
#     --mem 64G \
#     --job-name tensorqtl \
#     "${code_directory}/run_tensorqtl.sh" \
#         --genotype "${genotype}" \
#         --phenotype_list ${phenotype_list} \
#         --interactions "/oak/stanford/groups/smontgom/dnachun/projects/topmed/tensorqtl_inputs/Ready_for_running_inputs/interaction_terms/Eigen_Proteins_18modules_soft_thres_7_deep_split_3_me_0.25_min_size_20_4_removed_pcs.txt" \
#         --output_directory ${output_directory}
