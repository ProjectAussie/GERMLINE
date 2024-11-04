#!/usr/bin/env bash
set -euo pipefail

source ./testing_constants_and_functions.sh

output_file_prefix=output/haploid_test
germline_params_file="${output_file_prefix}_params.txt"
echo "Preparing germline params file ${germline_params_file}"
{
  echo 1
  echo test_dog.map
  echo single_dog_one_homoz_tract.ped
  echo ${output_file_prefix}
} > ${germline_params_file}

expected_output_file=one_dog_haploid_test_expected_without_snp_count.match

for bits in "${BITS_TO_TEST[@]}"; do
  log_file="${output_file_prefix}_${bits}_haploid_germline_log.txt"
  germline -haploid -silent -min_m 0.05 -err_hom 0 -err_het 0 -bits "${bits}" -w_extend \
    < "${germline_params_file}" \
    &> "${log_file}"
  compare_germline_output "${output_file_prefix}.match" "$expected_output_file" "$bits"
done
