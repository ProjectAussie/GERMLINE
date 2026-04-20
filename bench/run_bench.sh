#!/usr/bin/env bash
# Unattended peak-RSS benchmark for GERMLINE.
#
# Two modes:
#   BENCH_NATIVE=1 (default on Darwin) — build + run on the host (fast, OS-native sort).
#   BENCH_NATIVE=0                    — run inside the linux/amd64 builder container:
#     docker run --rm --platform=linux/amd64 -v "$PWD:/build" -w /build \
#       -e BENCH_NATIVE=0 --entrypoint bash germline-builder:latest -c "bash /build/bench/run_bench.sh"
#
# Inputs are supplied via env vars so no proprietary filenames are baked into the
# script. Defaults match the unchecked-in files in the repo root; override on any
# workstation that uses different names.
#
# Required / tunable env vars:
#   BENCH_LABEL       output-file prefix (default: bench)
#   BENCH_MAP         path to .map file                 (default: all_dogs_chr38.map)
#   BENCH_PED         path to .ped file                 (default: all_dogs_chr38.tail.ped)
#   BENCH_NEW         path to -new_samples list         (default: new_dogs.plinky)
#   BENCH_OLD         path to -samples_to_compare_to    (default: old_dogs.plinky)
#   BENCH_CHROM       chromosome id                     (default: 38)
#   BENCH_BITS        -bits value                       (default: 45)
#   BENCH_MIN_M       -min_m value                      (default: 0.5)
#   BENCH_REBUILD     1=make clean germline; 0=skip     (default: 1)
#   BENCH_NATIVE      1=native build, 0=in-container    (default: 1 on Darwin, 0 elsewhere)

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# default to native on macOS, in-container elsewhere
if [ -z "${BENCH_NATIVE:-}" ]; then
    if [ "$(uname -s)" = "Darwin" ]; then BENCH_NATIVE=1; else BENCH_NATIVE=0; fi
fi

if [ "${BENCH_NATIVE}" = "0" ]; then
    cd /build
else
    cd "${PROJECT_ROOT}"
fi

LABEL=${BENCH_LABEL:-bench}
MAP=${BENCH_MAP:-all_dogs_chr38.map}
PED=${BENCH_PED:-all_dogs_chr38.tail.ped}
NEW=${BENCH_NEW:-new_dogs.plinky}
OLD=${BENCH_OLD:-old_dogs.plinky}
CHROM=${BENCH_CHROM:-38}
BITS=${BENCH_BITS:-45}
MIN_M=${BENCH_MIN_M:-0.5}
REBUILD=${BENCH_REBUILD:-1}

# gnu time: gtime on macOS (brew coreutils), /usr/bin/time on Linux
if [ "${BENCH_NATIVE}" = "1" ] && [ "$(uname -s)" = "Darwin" ]; then
    TIME_BIN="${TIME_BIN:-gtime}"
    STAT_FMT="-f%z"   # macOS gstat syntax; we'll use wc -c instead for portability
    USE_GSTAT=1
else
    TIME_BIN="${TIME_BIN:-/usr/bin/time}"
    USE_GSTAT=0
fi

for f in "${MAP}" "${PED}" "${NEW}" "${OLD}"; do
    [ -f "${f}" ] || { echo "missing input: ${f}" >&2; exit 2; }
done

mkdir -p output

if [ "${REBUILD}" = "1" ]; then
    echo "=== build ==="
    make clean germline >/dev/null
fi

OUT_PREFIX="output/${LABEL}-chr${CHROM}"
TIME_LOG="${OUT_PREFIX}-time.log"
STDOUT_LOG="${OUT_PREFIX}-stdout.log"

rm -f "${OUT_PREFIX}".match "${OUT_PREFIX}".log "${TIME_LOG}" "${STDOUT_LOG}"

echo "=== run (${PED}, -bits ${BITS} -min_m ${MIN_M}, native=${BENCH_NATIVE}) ==="
"${TIME_BIN}" -v ./bin/germline \
    -chromosome "${CHROM}" \
    -haploid \
    -min_m "${MIN_M}" \
    -err_hom 0 \
    -err_het 0 \
    -bits "${BITS}" \
    -w_extend \
    -new_samples "${NEW}" \
    -samples_to_compare_to "${OLD}" <<EOF >"${STDOUT_LOG}" 2>"${TIME_LOG}"
1
${MAP}
${PED}
${OUT_PREFIX}
EOF

RSS_KB=$(awk '/Maximum resident set size/ {print $NF}' "${TIME_LOG}")
WALL=$(awk -F': ' '/Elapsed \(wall clock\)/ {print $NF}' "${TIME_LOG}")
if [ -f "${OUT_PREFIX}.match" ]; then
    MATCH_BYTES=$(wc -c < "${OUT_PREFIX}.match" | tr -d ' ')
    MATCH_LINES=$(wc -l < "${OUT_PREFIX}.match" | tr -d ' ')
else
    MATCH_BYTES=0
    MATCH_LINES=0
fi

echo "=== result (${LABEL}) ==="
printf "peak_rss_kb=%s\n" "${RSS_KB}"
printf "peak_rss_mb=%.1f\n" "$(awk -v x=${RSS_KB} 'BEGIN{print x/1024}')"
printf "wall_time=%s\n" "${WALL}"
printf "match_file_bytes=%s\n" "${MATCH_BYTES}"
printf "match_file_lines=%s\n" "${MATCH_LINES}"
