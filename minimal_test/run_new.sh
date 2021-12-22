#!/bin/bash

set -e


if [ "$1" == "old-new" ]; then
    ../germline -haploid -min_m 0.5 -err_hom 0 -err_het 0 -bits 90 -w_extend \
        -samples_to_compare_to 0000_old_dog_proxy_keys_germline_output.plinky \
        -new_samples 0000_new_dog_proxy_keys_germline_output.plinky \
        < germline_params_old_new
 #       -chromosome 1 \
#        -individual_outputs test_a \
fi

if [ "$1" == "new-new" ]; then
    ../germline -haploid -min_m 0.5 -err_hom 0 -err_het 0 -bits 90 -w_extend \
        -samples_to_compare_to 0000_new_dog_proxy_keys_germline_output.plinky \
        -new_samples 0000_new_dog_proxy_keys_germline_output.plinky \
        -chromosome 1 \
        -individual_outputs test_b \
        < germline_params_new_new
fi
