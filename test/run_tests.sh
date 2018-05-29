#!/bin/bash
set -euo pipefail
cd ..
make
cd test

bash test_single_dog_single_homoz_tract_w_extend.sh
bash test_two_dogs_single_homoz_tract_each_w_extend.sh
