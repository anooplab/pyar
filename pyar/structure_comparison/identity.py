"""Chemical identity provider interfaces and legacy adapters."""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from pyar.backends import babel
from pyar.structure_comparison.models import IdentityResult


@runtime_checkable
class ChemicalIdentityProvider(Protocol):
    """Service that identifies chemical species and compares their identity."""

    def identify(self, structure) -> IdentityResult:
        """Return a canonical chemical identity for a structure."""

    def same_identity(self, first, second) -> bool:
        """Return whether two structures have the same chemical identity."""


class OpenBabelIdentityProvider:
    """Adapter for PyAR's established OpenBabel/InChI identity decisions."""

    method = "openbabel-inchi"

    def identify(self, structure) -> IdentityResult:
        inchi = babel.make_inchi_string_from_xyz(structure)
        smiles = babel.make_smile_string_from_xyz(structure)
        if not inchi or not smiles:
            raise ValueError(
                f"Could not determine complete product identity from {structure}"
            )
        return IdentityResult(
            canonical_key=inchi,
            inchi=inchi,
            smiles=smiles,
            method=self.method,
        )

    def same_identity(self, first, second) -> bool:
        return self._coerce(first).canonical_key == self._coerce(second).canonical_key

    def _coerce(self, value) -> IdentityResult:
        if isinstance(value, IdentityResult):
            return value
        if isinstance(value, dict) and value.get("inchi"):
            return IdentityResult(
                canonical_key=value["inchi"],
                inchi=value["inchi"],
                smiles=value.get("smiles"),
                method=self.method,
            )
        return self.identify(value)
