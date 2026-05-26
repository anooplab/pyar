"""Shared helpers for xTB-backed backends."""


def build_xtb_command(executable, start_xyz_file, qc_params, opt_threshold=None):
    """Build an xTB command line from PyAR QC settings."""
    command = [executable, start_xyz_file]
    command.extend(xtb_parallel_args(qc_params))
    if opt_threshold is not None:
        command.extend(["-opt", opt_threshold])
    if qc_params.get("charge", 0) != 0:
        command.extend(["-chrg", str(qc_params["charge"])])
    if qc_params.get("multiplicity", 1) != 1:
        command.extend(["-uhf", str(qc_params["multiplicity"])])
    scftype = qc_params.get("scftype", "rhf")
    if qc_params.get("multiplicity", 1) == 1 and scftype != "rhf":
        command.append(f"-{scftype}")
    return command


def xtb_parallel_args(qc_params):
    """Return command-line arguments for xTB parallel execution."""
    nprocs = qc_params.get("nprocs")
    if nprocs is None:
        return []

    try:
        nprocs = int(nprocs)
    except (TypeError, ValueError):
        return []

    if nprocs < 1:
        return []

    return ["--parallel", str(nprocs)]
