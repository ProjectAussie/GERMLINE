#!/bin/bash

set -e

for f in `cat new_dogs.txt`;
do 
    mkdir -p test_rust_out_old_new/$f 
    mkdir -p test_rust_homoz_old_new/$f 
    mkdir -p test_rust_out_new_new/$f 
    mkdir -p test_rust_homoz_new_new/$f
done

./og_germline -haploid -min_m 0.5 -err_hom 0 -err_het 0 -bits 90 -w_extend \
    -samples_to_compare_to 0000_old_dog_proxy_keys_germline_output.plinky \
    -new_samples 0000_new_dog_proxy_keys_germline_output.plinky \
    < germline_params_old_new

./split_germline_output_for_new_dog.bin \
    --chromosome 1 \
    --new-dog-list new_dogs.txt \
    --pair-output-directory test_rust_out_old_new \
    --homoz-output-directory test_rust_homoz_old_new \
    < test_out_old_new.match

./og_germline -haploid -min_m 0.5 -err_hom 0 -err_het 0 -bits 90 -w_extend \
    -samples_to_compare_to 0000_new_dog_proxy_keys_germline_output.plinky \
    -new_samples 0000_new_dog_proxy_keys_germline_output.plinky \
    < germline_params_new_new

./split_germline_output_for_new_dog.bin \
    --chromosome 1 \
    --new-dog-list new_dogs.txt \
    --include-new-new-dog-pairs --include-new-dog-homoz \
    --pair-output-directory test_rust_out_new_new \
    --homoz-output-directory test_rust_homoz_new_new \
    < test_out_new_new.match
