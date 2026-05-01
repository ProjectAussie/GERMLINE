#!/usr/bin/env bash
# SCICO-1241: -individual_outputs in haploid mode opens one ofstream per
# (single_id) — not per Individual — so the userspace buffer for a per-dog
# match/homoz file can never interleave with another buffer pointed at the
# same path. This test exercises the path and asserts every record in every
# per-dog file is well-formed (NF==7 for match files, NF==3 for homoz files).
#
# A small dataset like this won't necessarily flush the 8 KiB filebuf often
# enough to land on the cut-mid-newline boundary that produced the prod
# corruption — that lives in the prod-scale benchmark on the dev server.
# This test is here to (a) keep -individual_outputs covered by CI at all
# (it had no coverage before SCICO-1241) and (b) catch obvious schema
# regressions in the per-dog output writers in Match.cpp.

set -euo pipefail

source ./testing_constants_and_functions.sh

output_dir=output/individual_outputs_haploid
rm -rf "$output_dir"
mkdir -p "$output_dir"

output_file_prefix=${output_dir}/run
germline_params_file="${output_file_prefix}_params.txt"
{
  echo 1
  echo test_dog.map
  echo individual_outputs_haploid.ped
  echo "${output_file_prefix}"
} > "${germline_params_file}"

log_file="${output_file_prefix}_germline_log.txt"
germline -haploid -silent -bits 21 -min_m 0.05 -err_hom 0 -err_het 0 -w_extend \
  -new_samples individual_outputs_new_samples.txt \
  -samples_to_compare_to individual_outputs_old_samples.txt \
  -individual_outputs "${output_dir}" \
  -chromosome 1 \
  < "${germline_params_file}" \
  &> "${log_file}"

match_dir="${output_dir}/dog_level_match_files"
homoz_dir="${output_dir}/dog_level_homoz_files"

if [ ! -d "${match_dir}" ]; then
  print_red "missing per-dog match dir ${match_dir}"
  exit 1
fi

# Every new dog must have at least one per-dog file (chr01.tsv).
for sid in 3201001 3201002; do
  for sub in dog_level_match_files dog_level_homoz_files; do
    f="${output_dir}/${sub}/${sid}/chr01.tsv"
    if [ ! -f "$f" ]; then
      print_red "missing expected output file ${f}"
      exit 1
    fi
  done
done

# Schema invariant: match TSVs are 7 fields, homoz TSVs are 3 fields. No
# empty rows. The bug that landed in v1.7 produced both NF=13 (two records
# concatenated mid-buffer-flush) and NF=0 (orphaned trailing newlines).
fail=0
while IFS= read -r f; do
  bad=$(awk 'NF != 7 { print NR ":" NF; exit }' "$f")
  if [ -n "$bad" ]; then
    print_red "match file ${f} has malformed row (line:nf = ${bad})"
    fail=1
  fi
done < <(find "${match_dir}" -type f -name '*.tsv')

while IFS= read -r f; do
  bad=$(awk 'NF != 3 { print NR ":" NF; exit }' "$f")
  if [ -n "$bad" ]; then
    print_red "homoz file ${f} has malformed row (line:nf = ${bad})"
    fail=1
  fi
done < <(find "${homoz_dir}" -type f -name '*.tsv')

if [ "$fail" -ne 0 ]; then
  exit 1
fi

print_green "individual_outputs (haploid): all per-dog rows well-formed"
