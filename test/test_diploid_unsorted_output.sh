#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
source ./testing_constants_and_functions.sh

echo "=== Diploid + w_extend + -unsorted_output (opt-in nondeterministic) tests ==="

# With -unsorted_output the .match file is written in hash-map iteration
# order, which is not stable across runs. The file still contains the same
# set of records, so comparing after a byte-level sort on each side is the
# correct equivalence check for this path.

output_prefix=output/diploid_w_extend_unsorted
params_file="${output_prefix}_params.txt"
{
  echo 1
  echo CEU.22.map
  echo CEU.22.ped
  echo ${output_prefix}
} > ${params_file}

germline -silent -unsorted_output -bits 50 -min_m 1 -err_hom 2 -err_het 0 -w_extend \
  < "${params_file}" &> "${output_prefix}_log.txt"

if diff -q <(LC_ALL=C sort expected_w_extend.match) <(LC_ALL=C sort "${output_prefix}.match"); then
  print_green "Diploid + w_extend + -unsorted_output (unrestricted): PASS"
else
  print_red "Diploid + w_extend + -unsorted_output (unrestricted): FAIL"
  exit 1
fi

output_prefix=output/restricted_w_extend_unsorted
params_file="${output_prefix}_params.txt"
{
  echo 1
  echo CEU.22.map
  echo CEU.22.ped
  echo ${output_prefix}
} > ${params_file}

germline -silent -unsorted_output -bits 50 -min_m 1 -err_hom 2 -err_het 0 -w_extend \
  -samples_to_compare_to old_humans -new_samples new_humans \
  < "${params_file}" &> "${output_prefix}_log.txt"

if diff -q <(LC_ALL=C sort restricted_w_extend.match) <(LC_ALL=C sort "${output_prefix}.match"); then
  print_green "Diploid + w_extend + -unsorted_output (restricted old vs new): PASS"
else
  print_red "Diploid + w_extend + -unsorted_output (restricted old vs new): FAIL"
  exit 1
fi
