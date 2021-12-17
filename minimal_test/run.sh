./og_germline -haploid -min_m 0.5 -err_hom 0 -err_het 0 -bits 90 -w_extend \
    -samples_to_compare_to 0000_old_dog_proxy_keys_germline_output.plinky \
    -new_samples 0000_new_dog_proxy_keys_germline_output.plinky \
    < germline_params_old_new

./og_germline -haploid -min_m 0.5 -err_hom 0 -err_het 0 -bits 90 -w_extend \
    -samples_to_compare_to 0000_new_dog_proxy_keys_germline_output.plinky \
    -new_samples 0000_new_dog_proxy_keys_germline_output.plinky \
    < germline_params_new_new
