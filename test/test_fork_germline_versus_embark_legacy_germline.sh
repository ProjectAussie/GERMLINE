#!/bin/bash

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
  for bits in 1 41; do
    echo "BITS: ${bits}"
    ${germline_bin} -haploid -min_m 0.5 -err_hom 0 -err_het 0 -bits ${bits} -w_extend < ${germline_params_file} 2>&1 | tee haploid_germline_log.txt || echo "Germline exit code: $?"
    awk '{$2=$2};1' haploid_test.match | cut -f5-7 -d' ' > haploid_test.match.bed
    awk -v bits=${bits} '{ print bits, $1, $2, $3 }' haploid_test.match.bed >> ${output_file}
  done
}

# The freshly compiled version of germline is one directory up, where we ran `make`
echo 'Testing embark fork of germline with -haploid -w_extend on a range of bits values'
echo 'bits chr start end' > compiled_germline_haploid_match_tracts_by_bits_parameter.txt
report_matches_for_bits_haploid ../germline compiled_germline_haploid_match_tracts_by_bits_parameter.txt

# The embark distribution of germline is in the PATH on an embark EC2
# echo 'Testing embark legacy germline (already in the PATH on an embark EC2) with -haploid -w_extend on a range of bits values'
# echo 'bits chr start end' > embark_germline_haploid_match_tracts_by_bits_parameter.txt
# report_matches_for_bits_haploid germline embark_germline_haploid_match_tracts_by_bits_parameter.txt

report_matches_for_bits_homoz_only () {
  germline_bin=$1
  output_file=$2

  germline_params_file=homoz_test_params.txt
  echo "Preparing germline params file ${germline_params_file}"
  echo 1 > ${germline_params_file}
  echo test_dog.map >> ${germline_params_file}
  echo test_dog.ped >> ${germline_params_file}
  echo homoz_test >> ${germline_params_file}

  for bits in 1 11 21 31 41 51 61 71 81 91 101 111 121; do
    ${germline_bin} -homoz-only -min_m 0.5 -err_hom 0 -err_het 0 -bits ${bits} -w_extend < ${germline_params_file} 2>&1 | tee homoz_germline_log.txt || echo "Germline exit code: $?"
    cat homoz_test.match >> ${output_file}
  done
}

# echo 'Testing embark fork of germline with -homoz-only -w_extend on a range of bits values'
# rm -f compiled_germline_homoz_match_tracts.match
# report_matches_for_bits_homoz_only ../germline compiled_germline_homoz_match_tracts.match

# echo 'Testing embark legacy germline (already in the PATH on an embark EC2) with -homoz-only -w_extend on a range of bits values'
# rm -f embark_germline_homoz_match_tracts.match
# report_matches_for_bits_homoz_only germline embark_germline_homoz_match_tracts.match

# echo 'DIFF OF FORK-GERMLINE VS. LEGACY GERMLINE IN HAPLOID MATCH TRACTS:'
# diff compiled_germline_haploid_match_tracts_by_bits_parameter.txt embark_germline_haploid_match_tracts_by_bits_parameter.txt

# echo 'DIFF OF FORK-GERMLINE VS. LEGACY GERMLINE IN HOMOZ MATCH TRACTS:'
# diff compiled_germline_homoz_match_tracts.match embark_germline_homoz_match_tracts.match
