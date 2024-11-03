#!/usr/bin/env bash

print_green () {
  tput setaf 2
  echo "$@"
  tput sgr0
}

print_red () {
  tput setaf 1
  echo "$@"
  tput sgr0
}

compare_germline_output() {
  local match_file=$1
  local expected_output_file=$2
  local bits=$3
  local actual_output_file="${match_file/.match/_bits_${bits}.match}"

  # The snp-count column (column 6) is inconsistent. This is a bug we'd like to fix in a future version. See the readme.
  # Sorting is important here because the order of match tracts coming out of germline is not consistent
  cut -f1-5,7-11 "$match_file" | sort > "$actual_output_file"
  if cmp "$actual_output_file" "$expected_output_file"; then
    print_green "bits=${bits} output file ($actual_output_file) matches expected outputs file $expected_output_file"
  else
    print_red "bits=${bits} output file ($actual_output_file) DOES NOT MATCH expected outputs file $expected_output_file"
  fi
}

# Our smaller match interval runs from SNP 6 to SNP 189, so we test -bits in increments of 20 up to 81.
# Germline finds the homoz match tract at -bits up to (including) 94. At bits 95, the first slice is SNPs 1-95, the second slice is SNPs 96-190, and the third slice is SNPs 191-200.
# Since none of these slices are fully homozygous, Germline doesn't find the homoz tract to extend from.
export BITS_TO_TEST=(11 21 41 61 81)
