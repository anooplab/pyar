"""Explicit atom-contact and operational-shell reporting tests."""

import numpy as np
import pytest

from pyar.core.molecule import Molecule
from pyar.microsolvation.contacts import assess_contacts
from pyar.microsolvation.placement import assess_shell
from pyar.microsolvation.request import MicrosolvationRequest
from pyar.persistence import read_json
from pyar.workflows.microsolvation import microsolvate


def methane():
    return Molecule(["C", "H", "H", "H", "H"],
                    [[0, 0, 0], [.629, .629, .629], [-.629, -.629, .629],
                     [-.629, .629, -.629], [.629, -.629, -.629]],
                    name="methane", title="methane", charge=0, multiplicity=1)


def water():
    return Molecule(["O", "H", "H"], [[0, 0, 0], [.9572, 0, 0], [-.239, .927, 0]],
                    name="water", title="water", charge=0, multiplicity=1)


def request(solute=None, groups=(("carbon", (0,)),), count=2):
    return MicrosolvationRequest(solute or methane(), water(), count, 4, {}, 2, None,
                                 32, 1.4, 3.5, site_groups=groups, contact_distance=3.5)


def candidate():
    req = request()
    molecule = req.solute.merged_with(req.solvent).merged_with(req.solvent)
    molecule.coordinates = np.vstack((req.solute.coordinates,
        req.solvent.coordinates + [3, 0, 0], req.solvent.coordinates + [4, 3, 0]))
    molecule.solute_atom_indices = tuple(range(5))
    molecule.solvent_fragments = ((5, 6, 7), (8, 9, 10))
    return req, molecule


def test_contact_evidence_counts_unique_molecules_and_overlapping_named_sites():
    req, molecule = candidate()
    req = request(groups=(("carbon", (0,)), ("overlap", (0,))))
    result = assess_contacts(molecule, req)
    assert result["status"] == "measured"
    assert [item["name"] for item in result["sites"]] == ["carbon", "overlap"]
    assert result["sites"][0]["occupancy_count"] == 1
    assert result["sites"][1]["occupancy_count"] == 1
    assert len(result["sites"][1]["occupied_by_solvents"]) == len(
        set(result["sites"][1]["occupied_by_solvents"]))
    assert result["whole_solute_first_shell_count"] == 1
    assert result["whole_solute_second_shell_count"] == 1
    assert result["solvent_molecules"][1]["shell"] == "second"
    assert result["interpretation"].startswith("geometric contacts")
    assert result["solute_connectivity"]["status"] in {"unchanged", "changed"}


def test_contact_model_tracks_post_relaxation_migration():
    req, molecule = candidate()
    molecule.coordinates[5:] += [20, 0, 0]
    result = assess_contacts(molecule, req)
    assert result["whole_solute_first_shell_count"] == 0
    assert all(site["occupancy_count"] == 0 for site in result["sites"])
    assert all(item["shell"] == "outer" for item in result["solvent_molecules"])


def test_contact_model_reports_unknown_when_original_solute_mapping_changes():
    req, molecule = candidate()
    molecule.atoms_list[0] = "N"
    assert assess_contacts(molecule, req)["status"] == "unknown"


def test_contact_metrics_are_embedded_in_the_existing_coverage_record():
    req, molecule = candidate()
    coverage = assess_shell(molecule, req)
    assert coverage["contact_analysis"]["whole_solute_first_shell_count"] == 1
    assert coverage["contact_analysis"]["whole_solute_second_shell_count"] == 1
    assert coverage["coverage_fraction"] >= 0


def test_contact_metrics_survive_v2_step_commit_and_completed_rerun(tmp_path):
    req = request(count=1)
    first = microsolvate(req, output=tmp_path / "contacts")
    state = read_json(tmp_path / "contacts" / "state.json")
    selection = read_json(tmp_path / "contacts" / "steps" / "step_001" / "selection.json")
    stored = selection["selected"][0]["coverage"]["contact_analysis"]
    assert stored["status"] == "measured"
    assert state["coverage_history"][0][0]["contact_analysis"] == stored
    second = microsolvate(req, output=tmp_path / "contacts")
    assert first.selected_paths == second.selected_paths


def test_site_group_indices_and_names_are_validated():
    with pytest.raises(ValueError, match="distinct non-empty"):
        request(groups=(("duplicate", (0,)), ("duplicate", (1,))))
    with pytest.raises(ValueError, match="distinct non-empty"):
        request(groups=(("repeated-index", (0, 0)),))
    with pytest.raises(ValueError, match="out-of-range"):
        request(groups=(("bad", (5,)),))
    with pytest.raises(ValueError, match="contact reporting requires"):
        MicrosolvationRequest(methane(), water(), 1, site_groups=(("x", (0,)),))


def test_site_groups_are_canonicalized_for_state_and_request_ownership():
    req = request(groups=[("  proton donor  ", ["0", 1])])

    assert req.site_groups == (("proton donor", (0, 1)),)
    assert req.to_state_dict()["site_contact_policy"]["groups"] == [
        {"name": "proton donor", "solute_atom_indices": [0, 1]}
    ]


def test_site_group_validation_keeps_first_error_order():
    with pytest.raises(ValueError, match="distinct non-empty"):
        MicrosolvationRequest(
            methane(), water(), 1,
            site_groups=(("duplicate", (0,)), ("duplicate", (1,))),
            contact_distance=0.0,
        )
    with pytest.raises(ValueError, match="out-of-range"):
        MicrosolvationRequest(
            methane(), water(), 1,
            site_groups=(("outside", (5,)),),
            contact_distance=0.0,
        )
    with pytest.raises(ValueError, match="--contact-distance requires at least one"):
        MicrosolvationRequest(methane(), water(), 1, contact_distance=3.0)


def test_anionic_substrate_can_be_reported_without_assuming_donor_acceptor_roles():
    acetate = Molecule(["C", "C", "O", "O", "H", "H", "H"],
                       [[0, 0, 0], [1.5, 0, 0], [-.7, 1.1, 0], [-.7, -1.1, 0],
                        [0, 0, 1], [2, 0, 1], [2, 0, -1]],
                       name="acetate", title="acetate", charge=-1, multiplicity=1)
    req = request(solute=acetate, groups=(("carboxylate", (2, 3)),))
    assert req.to_state_dict()["site_contact_policy"]["groups"][0]["solute_atom_indices"] == [2, 3]
    molecule = acetate.merged_with(water())
    molecule.coordinates = np.vstack((acetate.coordinates, water().coordinates + [1.7, 1.1, 0]))
    molecule.solute_atom_indices = tuple(range(len(acetate)))
    molecule.solvent_fragments = (tuple(range(len(acetate), len(molecule))),)
    result = assess_contacts(molecule, req)
    assert result["sites"][0]["occupancy_count"] == 1
