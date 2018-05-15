#!/bin/bash
cd ..
make
cd test

source ./testing_constants_and_functions.sh
mkdir -p test_outputs

germline_params_file=test_outputs/haploid_test_params.txt
echo "Preparing germline params file ${germline_params_file}"
echo 1 > ${germline_params_file}
echo test_dog.map >> ${germline_params_file}
echo single_dog_one_homoz_tract.ped >> ${germline_params_file}
echo haploid_test >> ${germline_params_file}

for bits in "${BITS_TO_TEST[@]}"; do
  echo "BITS: ${bits}"
  ../germline -haploid -silent -min_m 0.05 -err_hom 0 -err_het 0 -bits ${bits} -w_extend < ${germline_params_file} 2>&1 | tee haploid_germline_log.txt || echo "Germline exit code: $?"
  cut -f1-5,7-11 haploid_test.match | sort > test_outputs/one_dog_haploid_test_bits_${bits}.match
  cmp test_outputs/one_dog_haploid_test_bits_${bits}.match one_dog_haploid_test_expected_without_snp_count.match && print_green "bits=${bits} output file (test_outputs/one_dog_haploid_test_bits_${bits}.match) matches expected outputs file one_dog_haploid_test_expected_without_snp_count.match" || (print_red "bits=${bits} output file (test_outputs/one_dog_haploid_test_bits_${bits}.match) DOES NOT MATCH expected outputs file one_dog_haploid_test_expected_without_snp_count.match" && exit 1)
done
