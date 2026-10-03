"""Repeat small real-backend scan smoke checks; these are not TS benchmarks.

Run from a PyAR installation with the requested external executable available.
For ORCA GFN2-xTB, configure ORCA's XTBEXE setting before invoking this script.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from pyar.workflows.scan_bond import run_scan_bond


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--backend", action="append", choices=("xtb-gfn2", "xtb-gxtb", "orca-gfn2"))
    args = parser.parse_args(argv)
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    fragment = root / "H.xyz"
    fragment.write_text("1\nH atom for H2 scan smoke test\nH 0 0 0\n")
    records = []
    for backend in args.backend or ["xtb-gfn2"]:
        software, model = backend.split("-")
        params = {"software": software, "nprocs": 1, "opt_cycles": 100}
        params.update({"method": "GFN2-xTB"} if software == "orca" else {"xtb_model": model})
        output = root / backend
        result = run_scan_bond(fragment, fragment, [0, 0], 1, params, output,
                               scan_points=3, scan_end=0.74, through="scan")
        orientation = result["results"][0]
        record = {"backend": backend, "status": result["status"],
                  "scan_status": orientation["scan_status"], "output": str(output)}
        if result["status"] == "complete":
            profile = json.loads(Path(orientation["scan_profile_json"]).read_text())["points"]
            energies = [point["energy_hartree"] for point in profile]
            record.update(points=len(profile), energies_hartree=energies,
                          finite_energies=bool(np.isfinite(energies).all()))
            if len(profile) != 3 or not record["finite_energies"]:
                raise ValueError(f"Invalid scan result: {record}")
        else:
            record["error"] = orientation.get("error")
        records.append(record)
    (root / "smoke_results.json").write_text(json.dumps(records, indent=2))
    print(json.dumps(records, indent=2))
    return int(any(record["status"] != "complete" for record in records))


if __name__ == "__main__":
    raise SystemExit(main())
