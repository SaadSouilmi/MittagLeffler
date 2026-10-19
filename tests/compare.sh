#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
Stage=${1:-all}
Results=results/real-arithmetic
Baseline=$Results/baseline
CPU=${ML_CPU:-0}
mkdir -p "$Results"
if [[ "$Stage" == build || "$Stage" == all ]]; then
	if [[ ! -d "$Baseline/src" ]]; then
		mkdir -p "$Baseline"
		git archive a17defe | tar -x -C "$Baseline"
	fi
	git rev-parse a17defe > "$Results/baseline-commit.txt"
	lscpu > "$Results/cpu.txt"
	"${CXX:-g++}" --version > "$Results/compiler.txt"
	printf '%s\n' "$CPU" > "$Results/cpu-affinity.txt"
	for Variant in baseline real; do
		Root=.
		if [[ "$Variant" == baseline ]]; then Root=$Baseline; fi
		Flags=(-std=gnu++20 -O3 -march=native -Wall -Wextra -Wno-unused-function -Wno-missing-field-initializers -pipe)
		if [[ "$Variant" == baseline ]]; then Flags+=(-include complex.h); fi
		"${CXX:-g++}" "${Flags[@]}" "$Root/src/exp/throughput.cpp" -o "$Results/$Variant-throughput" -lm
		for Test in correctness initialization; do
			"${CXX:-g++}" "${Flags[@]}" -I "$Root/src" "tests/$Test.cpp" -o "$Results/$Variant-$Test" -lm
		done
		objdump -d --no-show-raw-insn -C "$Results/$Variant-initialization" > "$Results/$Variant-initialization.asm"
		objdump -d --no-show-raw-insn -C "$Results/$Variant-throughput" > "$Results/$Variant-throughput.asm"
	done
fi
if [[ "$Stage" == run || "$Stage" == all ]]; then
	for Variant in baseline real; do
		"$Results/$Variant-correctness" > "$Results/$Variant-samples.csv" 2> "$Results/$Variant-correctness.log"
	done
	for Run in 1 2 3 4 5; do
		Variants=(baseline real)
		if (( Run % 2 == 0 )); then Variants=(real baseline); fi
		for Variant in "${Variants[@]}"; do
			printf 'Run %s: %s\n' "$Run" "$Variant"
			taskset -c "$CPU" "$Results/$Variant-throughput" > "$Results/$Variant-throughput-$Run.txt"
			taskset -c "$CPU" "$Results/$Variant-initialization" > "$Results/$Variant-initialization-$Run.csv"
		done
	done
fi
if [[ "$Stage" == analyze || "$Stage" == all ]]; then
	uv run --python 3.12 tests/compare.py "$Results"
fi
