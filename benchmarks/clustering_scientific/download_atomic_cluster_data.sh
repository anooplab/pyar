#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ccd_dir="$root/data/cambridge_cluster_database"
qcd_dir="$root/data/qcd_au13"
mkdir -p "$ccd_dir" "$qcd_dir"

curl -fL --retry 3 \
  'https://www-wales.ch.cam.ac.uk/~wales/CCD/LJ13.tar.bz2' \
  -o "$ccd_dir/LJ13.tar.bz2"

# The QCD publication is archived in NOMAD. Query its Au13 entries from the
# published dataset and request only geometry, total energy, and provenance.
curl -fLsS --retry 3 --max-time 180 \
  -H 'Content-Type: application/json' -X POST \
  'https://nomad-lab.eu/prod/v1/api/v1/entries/archive/download/query' \
  --data-binary '{"query":{"datasets.dataset_id":"vDrK65kaR1myz2YRXjyZ3A","results.material.chemical_formula_hill":"Au13"},"required":{"metadata":{"entry_id":"*","mainfile":"*","upload_id":"*","license":"*"},"results":{"properties":{"geometry_optimization":{"final_force_maximum":"*","final_energy_difference":"*","final_displacement_maximum":"*"}}},"run":{"system[-1]":{"atoms":"*"},"calculation[-1]":{"energy":{"total":"*"}}}}}' \
  -o "$qcd_dir/nomad_archive.zip"

echo 'Downloaded source data. Verify checksums against data/SOURCES.md.'
