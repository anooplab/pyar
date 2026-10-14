# Legacy grow fixture

`grow_interrupted/` preserves raw snapshot 001 from PR 0 at baseline
8ab5bcca43d43fbcf69a76644fad79e1c233b37e. It is the committed first step
of `test_interrupted_restart_reuses_pool`: H2 seed + H2 monomer, three additions,
mocked addition engine, no invented energies. Request, state and selected XYZ
bytes are unmodified; all references are relative, with no machine paths.
The fixture can be copied and loaded directly by the original v1 owner.
Source file hashes are in docs/run_layout_v2/baseline_inventory.json.
This is storage/scientific-input regression evidence, not a chemical benchmark.
