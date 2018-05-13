#!/bin/bash
cd ..
make
cd test

print_green () {
  tput setaf 2
  echo $@
  tput setaf 7
}

print_red () {
  tput setaf 1
  echo $@
  tput setaf 7
}

report_matches_for_bits_haploid () {
  germline_bin=$1
  output_file=$2

  germline_params_file=haploid_test_params.txt
  echo "Preparing germline params file ${germline_params_file}"
  echo 1 > ${germline_params_file}
  echo test_dog.map >> ${germline_params_file}
  echo test_dog.ped >> ${germline_params_file}
  echo haploid_test >> ${germline_params_file}

  # for bits in 1 11 21 31 41 51 61 71 81 91 101 111 121; do
  for bits in 1 21 41 61 81; do
    echo "BITS: ${bits}"
    ${germline_bin} -haploid -silent -min_m 0.5 -err_hom 0 -err_het 0 -bits ${bits} -w_extend < ${germline_params_file} 2>&1 | tee haploid_germline_log.txt || echo "Germline exit code: $?"
    awk '{$2=$2};1' haploid_test.match | cut -f5-7 -d' ' > haploid_test.match.bed
    awk -v bits=${bits} '{ print bits, $1, $2, $3 }' haploid_test.match.bed >> ${output_file}
  done
}

# The freshly compiled version of germline is one directory up, where we ran `make`
echo 'Testing embark fork of germline with -haploid -w_extend on a range of bits values'
echo 'bits chr start end' > compiled_germline_haploid_match_tracts_by_bits_parameter.txt
report_matches_for_bits_haploid ../germline compiled_germline_haploid_match_tracts_by_bits_parameter.txt

cat compiled_germline_haploid_match_tracts_by_bits_parameter.txt | tail -n +2 | cut -d' ' -f2-4 | uniq > unique_haploid_outputs.txt
unique_lines_in_haploid_output=$(wc -l < unique_haploid_outputs.txt | xargs)
if [[ "${unique_lines_in_haploid_output}" = 1 ]]; then
  tput setaf 2
  echo "-haploid -w_extend test passes. All germline -haploid -w_extend runs produced the same match interval."
  cat compiled_germline_haploid_match_tracts_by_bits_parameter.txt
  tput setaf 7
else
  tput setaf 1
  echo "-haploid -w_extend test fails. Expected all germline -haploid -w_extend runs to produce the same start and end coordinates. Got:"
  cat compiled_germline_haploid_match_tracts_by_bits_parameter.txt
  tput setaf 7
  exit 1
fi

report_matches_for_bits_homoz_only () {
  germline_bin=$1
  output_file=$2

  germline_params_file=homoz_test_params.txt
  echo "Preparing germline params file ${germline_params_file}"
  echo 1 > ${germline_params_file}
  echo test_dog.map >> ${germline_params_file}
  echo test_dog.ped >> ${germline_params_file}
  echo homoz_test >> ${germline_params_file}

  for bits in 1 21 41 61 81; do
    ${germline_bin} -homoz-only -silent -min_m 0.5 -err_hom 0 -err_het 0 -bits ${bits} -w_extend < ${germline_params_file} 2>&1 | tee homoz_germline_log.txt || echo "Germline exit code: $?"
    cat homoz_test.match >> ${output_file}
  done
}

echo 'Testing embark fork of germline with -homoz-only -w_extend on a range of bits values'
rm -f compiled_germline_homoz_match_tracts.match
report_matches_for_bits_homoz_only ../germline compiled_germline_homoz_match_tracts.match

cat compiled_germline_homoz_match_tracts.match | uniq > unique_homoz_outputs.txt
unique_lines_in_homoz_output=$(wc -l < unique_homoz_outputs.txt | xargs)
if [[ "${unique_lines_in_haploid_output}" = 1 ]]; then
  tput setaf 2
  echo "-homoz-only -w_extend test passes. All germline -homoz-only -w_extend runs produced the same match interval."
  cat compiled_germline_homoz_match_tracts.match
  tput setaf 7
else
  tput setaf 1
  echo "-homoz-only -w_extend test fails. Expected all germline -homoz-only -w_extend runs to produce the same start and end coordinates. Got:"
  cat compiled_germline_homoz_match_tracts.match
  tput setaf 7
  exit 1
fi

