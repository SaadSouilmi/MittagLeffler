#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
./build.sh
"${CXX:-g++}" -std=gnu++20 -O3 -march=native -Wall -Wextra \
	-Wno-unused-function -Wno-missing-field-initializers -pipe \
	src/exp/accuracy.cpp -o bin/accuracy -lm
printf 'ML_METADATA_BEGIN\n'
lscpu | sed -n '/Model name/p'
"${CXX:-g++}" --version | head -1
printf 'CPU_AFFINITY=0\nML_METADATA_END\n'
for Run in 1 2 3 4 5; do
	printf 'ML_BENCHMARK_RUN=%s\n' "$Run"
	taskset -c 0 ./bin/throughput
done
printf 'ML_SAMPLES_BEGIN\n'
./bin/accuracy
printf 'ML_SAMPLES_END\n'
