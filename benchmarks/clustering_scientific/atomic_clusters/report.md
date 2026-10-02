# Atomic cluster benchmark results

Run artifacts were generated with the repository script against 1,510 CCD
LJ13 minima and 30 QCD/NOMAD Au13 records. The complete 30-condition results,
individual input structures, source IDs, force diagnostics, and checksums are
saved under
[`../runs/atomic_clusters/`](../runs/atomic_clusters/).

## Outcome

The table orders conditions by the benchmark's normalized mean nearest
pair-spectrum RMS (lower means the selected set lies closer to more of the
input structures under this descriptor proxy):

| Dataset | Requested algorithm | Feature | Selected | Mean normalized pair-spectrum RMS | Lowest-energy entry retained |
|---|---|---:|---:|---:|---|
| LJ13 | DBSCAN | MBTR | 12 | 0.05097 | Yes |
| LJ13 | OPTICS | SOAP | 12 | 0.05216 | Yes |
| LJ13 | auto | MBTR | 12 | 0.05228 | Yes |
| Au13 | auto | SOAP | 12 | 0.02050 | Yes |
| Au13 | agglomerative | distance histogram | 12 | 0.02060 | Yes |
| Au13 | agglomerative | SOAP | 12 | 0.02075 | Yes |

Every tested condition retained the lowest-energy source structure under the
energies supplied in that dataset. The best feature/algorithm combination
differs between LJ13 and Au13, which argues against selecting a universal
atomic-cluster policy from these two pilots. These are descriptive scores, not
verified recovery rates: sorted pair-distance spectra can collide for
non-congruent geometries, and no independent basin labels are available.

For LJ13, auto/HDBSCAN with SOAP produced only four cluster minima, so the
12-seed budget could not be filled. Agglomerative and density-based methods
produced enough candidates in the other conditions. OPTICS with the distance
histogram completed but emitted a recorded scikit-learn divide-by-zero
`RuntimeWarning`; treat that result as requiring follow-up before relying on
that configuration.

QCD includes a final-force diagnostic for each Au13 entry. Converted values
range up to 0.0598 eV/angstrom. The entries are not filtered by one force
threshold or a single documented electronic-structure protocol in this
benchmark, so source-energy ranks are reported within this archive only.

The requested `maxmin` setting follows the implemented production sequence:
automatic clustering first, then max-min trimming only if cluster minima
exceed the seed budget. It therefore shares the `auto` clustering result; it
is not a standalone max-min partitioner.

## Reproduction

```bash
bash benchmarks/clustering_scientific/download_atomic_cluster_data.sh
python -m pyar.scripts.atomic_cluster_benchmark \
  --lj-archive benchmarks/clustering_scientific/data/cambridge_cluster_database/LJ13.tar.bz2 \
  --au13-archive benchmarks/clustering_scientific/data/qcd_au13/nomad_archive.zip \
  --output benchmarks/clustering_scientific/runs/atomic_clusters
```

Source descriptions, terms, unit conversion details, and SHA-256 values are
documented in [`../data/SOURCES.md`](../data/SOURCES.md). See
[`README.md`](README.md) for the benchmark's scope and limitations.
