#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p bin
"${CXX:-g++}" -std=gnu++20 -O3 -march=native -Wall -Wextra \
	-Wno-unused-function -Wno-missing-field-initializers -pipe \
	src/exp/throughput.cpp -o bin/throughput -lm
