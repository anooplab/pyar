# Atomic cluster benchmark sources

## Lennard-Jones 13

`cambridge_cluster_database/LJ13.tar.bz2` is the Cambridge Cluster Database
LJ13 archive from <https://www-wales.ch.cam.ac.uk/~wales/CCD/LJ13.tar.bz2>.
The CCD lists 1,510 minima and 29,007 transition states for this landscape;
the benchmark uses the minima (`min.data` and `extractedmin`). Energies are in
Lennard-Jones epsilon and coordinates in sigma. The CCD page does not state an
explicit reuse license for this archive, so this repository records the
source and checksum but does not assert a license. Consult the source terms
before redistributing the archive beyond this research repository.

SHA-256: `f915652c397a3c253985e6c88c3e21e965f6deaeb79c717627e1cf85895fd0ff`

## Gold 13

`qcd_au13/nomad_archive.zip` contains the 30 Au13 entries queried from the
Quantum Cluster Database's published NOMAD dataset, DOI
[`10.17172/NOMAD/2023.02.01-1`](https://doi.org/10.17172/NOMAD/2023.02.01-1),
using dataset ID `vDrK65kaR1myz2YRXjyZ3A`. The QCD paper is
[Müller et al., Scientific Data 10, 39 (2023)](https://www.nature.com/articles/s41597-023-02200-4).
The returned archive records declare **CC BY 4.0**. Attribute the QCD authors
and dataset, retain entry IDs and source mainfiles in `structures.csv`, and
link the license at <https://creativecommons.org/licenses/by/4.0/>.

Coordinates are NOMAD SI metres converted to angstrom; total energies are
joules converted to eV. The downloaded archive includes geometry-optimization
force/convergence metadata. It does not follow that all QCD entries are
independent minima under one identical computational protocol; they are
treated as source entries and never merged as ground-truth labels here.

SHA-256: `77e447eb139654542049a67106b2170ceed373b7f8d9dc63f5517be4b2c1dd83`

Re-download both sources with `bash benchmarks/clustering_scientific/download_atomic_cluster_data.sh`.

## Water hexamer minima

`water_clusters/W6_geoms_5.0_KCal-1hgztfv.txt` is the per-size W6 file linked
from the University of Washington's [Database of Water Cluster Minima](https://sites.uw.edu/wdbase/database-of-water-clusters/).
It contains 80 (H2O)6 minima within 5 kcal/mol of the reported putative
minimum, with energies and coordinates from the TTM2.1-F potential. The larger
database reports 4,948,959 minima for sizes n=3–30. This benchmark uses one
small fixed-size source subset and then deterministically selects 24 source
energy ranks spanning its range.

The source table does not state a separate reuse license for this file. The
upstream URL, file hash, and attribution are retained; check with the source
maintainers before redistributing these geometries outside this research
repository.

SHA-256: `c2f7f33ee5028c3029ed1271f110edd88d88dfacddc530bfce8e64c15483ae5e`

Re-download with `bash benchmarks/clustering_scientific/download_water_cluster_sample.sh`.
