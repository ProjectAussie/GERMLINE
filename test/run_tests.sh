#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT=$(realpath "$(dirname "$0")"/..)
export PATH=${REPO_ROOT}/bin:${PATH}
cd "$REPO_ROOT/test"
mkdir -p output

./test_single_dog_single_homoz_tract_w_extend.sh
./test_two_dogs_single_homoz_tract_each_w_extend.sh
