"""Compare saved one-direction water reference with production bidirectional rule."""
import json
import time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from ase.io import read
from pyar.structure_comparison.fragment_rmsd import FragmentRMSDComparator

HERE=Path(__file__).resolve().parent
source=HERE.parent / "runs/water_similarity"
frames=read(source / "structures.xyz", index=":")
molecules=[SimpleNamespace(atoms_list=a.get_chemical_symbols(), coordinates=a.positions) for a in frames]
reference=np.loadtxt(source / "reference_fragment_rmsd_upper_bound_angstrom.csv", delimiter=",")
left,right=np.triu_indices(len(frames),1)
order=np.argsort(abs(reference[left,right]-.75))[:3]
comparator=FragmentRMSDComparator(atom_mode="all",max_mappings=10000)
results=[]
for k in order:
    i,j=int(left[k]),int(right[k])
    started=time.monotonic()
    forward=comparator.compare(molecules[i],molecules[j])
    reverse=comparator.compare(molecules[j],molecules[i])
    row={"left":i,"right":j,"saved_distance":float(reference[i,j]),
         "forward":forward.distance,"reverse":reverse.distance,
         "production_distance":max(forward.distance,reverse.distance),
         "elapsed_seconds":time.monotonic()-started}
    results.append(row)
    Path(__file__).with_suffix(".json").write_text(json.dumps(results,indent=2)+"\n")
    print(row,flush=True)
