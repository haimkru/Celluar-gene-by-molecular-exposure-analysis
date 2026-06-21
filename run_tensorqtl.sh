#!/usr/bin/env bash

set -o xtrace -o nounset -o errexit

check_for_file() {
    argument_name="${1}"
    file_path="${2}"
    if [[ ${file_path} != "none" ]] && [[ ! -f ${file_path} ]]; then
        echo "Error: file ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

check_for_directory() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ${directory_path} != "none" ]] && [[ ! -d ${directory_path} ]]; then
        echo "Error: directory ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

options_array=(
    genotype
    phenotype_list
    interactions
    output_directory
)

eval "$(printf "%s\n" "${options_array[@]}" | xargs --replace=% echo "declare %=none;")"
longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g' | sed -e 's/$/:/')

# Parse command line arguments with getopt
arguments=$(getopt --options a --longoptions "${longoptions}" --name 'run_tensorqtl' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --genotype )
            genotype="${2}"; shift 2 ;;
        --phenotype_list )
            phenotype_list="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --interactions )
            interactions="${2}"; check_for_file "${1}" "${2}"; shift 2 ;;
        --output_directory )
            output_directory="${2}"; shift 2 ;;
        -- )
            shift; break;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

export PATH="${HOME}/.pixi/bin:/oak/stanford/groups/smontgom/dnachun/pixi/bin:${PATH}"

line_number=${SLURM_ARRAY_TASK_ID}
phenotype_bed="$(sed "${line_number}q; d" "${phenotype_list}")"
prefix=$(basename ${phenotype_bed%.bed})
covariates=${phenotype_bed//standardized_phenotypes/covariates}
covariates=${covariates//.bed/_covariates.txt}

if [[ ${prefix} =~ .*expression.* ]]; then
    window_size="1000000"
else
    window_size="5000"
    interactions=${interactions//Expression/Methylation}
fi

tensorqtl \
    "${genotype}" \
    "${phenotype_bed}" \
    ${prefix} \
    --covariates "${covariates}" \
    --interaction "${interactions}" \
    --mode cis_nominal \
    --maf_threshold 0.05 \
    --pval_threshold 1 \
    --seed 1 \
    --window ${window_size} \
    --output_dir ${output_directory}
