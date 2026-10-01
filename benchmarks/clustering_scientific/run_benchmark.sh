#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$root"
data="benchmarks/clustering_scientific/data/mpconf196gen"
runs="benchmarks/clustering_scientific/runs"
mkdir -p "$runs"

for system in CAMVES_I FGG55 WG01; do
    python -m pyar.scripts.scientific_clustering_benchmark \
        "$data/${system}_crest_conformers.xyz" \
        --output "$runs/$system" \
        --max-seeds 12 \
        --downstream-perturbation 0.10
done

python -m pyar.scripts.summarize_scientific_clustering \
    --runs "$runs/CAMVES_I" "$runs/FGG55" "$runs/WG01" \
    --output benchmarks/clustering_scientific/analysis
