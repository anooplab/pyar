"""Shared helpers for xTB-backed interfaces."""


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
