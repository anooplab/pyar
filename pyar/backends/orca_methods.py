"""Method keyword helpers shared by ORCA calculation paths."""

from __future__ import annotations


_ORCA_XTB_METHODS = {
    "gxtb": "g-xTB",
    "gfn0xtb": "GFN0-xTB",
    "xtb0": "GFN0-xTB",
    "gfnxtb": "GFN-xTB",
    "gfn1xtb": "GFN-xTB",
    "xtb1": "GFN-xTB",
    "gfn2xtb": "GFN2-xTB",
    "xtb2": "GFN2-xTB",
    "xtb": "GFN2-xTB",
    "gfnff": "GFN-FF",
    "xtbff": "GFN-FF",
}


def orca_method(method):
    """Return ``(ORCA method spelling, is_xTB)`` for a user method name."""
    supplied = str(method).strip()
    key = "".join(character for character in supplied.lower() if character.isalnum())
    canonical = _ORCA_XTB_METHODS.get(key)
    return (canonical, True) if canonical else (supplied, False)


def orca_method_keywords(qc_params, optimization_keyword=None):
    """Build the method part of an ORCA keyword line.

    ORCA's xTB methods do not take a Gaussian basis or the DFT-only RI,
    dispersion, and SCF accelerator keywords used by PyAR's DFT defaults.
    """
    method, is_xtb = orca_method(qc_params.get("method", "BP86"))
    if method == "g-xTB":
        pieces = ["!", "ExtOpt"]
        if optimization_keyword:
            pieces.append(optimization_keyword)
        return " ".join(pieces), True
    pieces = ["!", method]
    if not is_xtb and not qc_params.get("orca_builtin_method"):
        basis = qc_params.get("basis")
        if not basis:
            raise ValueError("ORCA DFT calculations require a basis set")
        pieces.extend([str(basis), "RI", "def2/J", "D3BJ", "KDIIS"])
    if optimization_keyword:
        pieces.append(optimization_keyword)
    return " ".join(pieces), is_xtb


def orca_external_method_block(qc_params):
    """Return the ORCA ``ProgExt`` block for external methods, if selected."""
    method, _ = orca_method(qc_params.get("method", "BP86"))
    if method != "g-xTB":
        return ""
    wrapper = qc_params.get("gxtb_wrapper")
    if not wrapper:
        raise ValueError("g-xTB with ORCA requires the path to the oet_gxtb wrapper")
    escaped_path = str(wrapper).replace('"', '\\"')
    return f'%method\n  ProgExt "{escaped_path}"\nend'
