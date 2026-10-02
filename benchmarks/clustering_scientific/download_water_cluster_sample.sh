#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data_dir="$root/data/water_clusters"
mkdir -p "$data_dir"
curl -fL --retry 3 \
  'https://sites.uw.edu/wdbase/files/2019/01/W6_geoms_5.0_KCal-1hgztfv.txt' \
  -o "$data_dir/W6_geoms_5.0_KCal-1hgztfv.txt"

echo 'Downloaded the W6 TTM2.1-F minima subset; verify against data/SOURCES.md.'
