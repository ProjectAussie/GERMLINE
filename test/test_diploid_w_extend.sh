#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
source ./testing_constants_and_functions.sh

echo "=== Diploid + w_extend tests ==="

# Test 1: Unrestricted diploid with w_extend
output_prefix=output/diploid_w_extend
params_file="${output_prefix}_params.txt"
{
  echo 1
  echo CEU.22.map
  echo CEU.22.ped
  echo ${output_prefix}
} > ${params_file}

germline -silent -bits 50 -min_m 1 -err_hom 2 -err_het 0 -w_extend \
  < "${params_file}" &> "${output_prefix}_log.txt"

sort expected_w_extend.match > output/expected_w_extend_sorted.match
sort "${output_prefix}.match" > output/diploid_w_extend_sorted.match
if diff -q output/expected_w_extend_sorted.match output/diploid_w_extend_sorted.match > /dev/null 2>&1; then
  print_green "Diploid + w_extend (unrestricted): PASS"
else
  print_red "Diploid + w_extend (unrestricted): FAIL"
  exit 1
fi

# Test 2: Restricted diploid with w_extend (old vs new)
output_prefix=output/restricted_w_extend
params_file="${output_prefix}_params.txt"
{
  echo 1
  echo CEU.22.map
  echo CEU.22.ped
  echo ${output_prefix}
} > ${params_file}

germline -silent -bits 50 -min_m 1 -err_hom 2 -err_het 0 -w_extend \
  -samples_to_compare_to old_humans -new_samples new_humans \
  < "${params_file}" &> "${output_prefix}_log.txt"

sort restricted_w_extend.match > output/restricted_w_extend_expected_sorted.match
sort "${output_prefix}.match" > output/restricted_w_extend_sorted.match
if diff -q output/restricted_w_extend_expected_sorted.match output/restricted_w_extend_sorted.match > /dev/null 2>&1; then
  print_green "Diploid + w_extend (restricted old vs new): PASS"
else
  print_red "Diploid + w_extend (restricted old vs new): FAIL"
  exit 1
fi
