from __future__ import annotations

import copy
import json

import pytest

from phenocycler import marker_taxonomy
from phenocycler.marker_registry import (
    DEFAULT_REGISTRY_PATH,
    load_registry,
    registry_from_dict,
)


def test_registry_covers_current_type_and_state_rosters():
    registry = load_registry()
    assert {marker.name for marker in registry.type_markers} == set(
        marker_taxonomy.TYPE
    )
    assert {marker.name for marker in registry.state_markers} == set(
        marker_taxonomy.PROCESS
    )
    assert len(registry.by_name) == len(registry.markers)
    assert len(registry.fingerprint) == 64


def test_registry_locks_broad_and_subtype_decisions():
    registry = load_registry()
    broad = {marker.name for marker in registry.broad_markers}
    assert len(broad) == 23
    assert "SST" not in broad
    assert registry.marker("SST").subtype_only
    assert registry.marker("SST").reference == "Pan_Cytokeratin"
    assert registry.marker("CD11b").reference == "E_cadherin"
    assert set(registry.reference_graph) == {
        marker.name for marker in registry.active_markers
    }


def test_familywise_allocations_are_fixed_at_one_percent_per_node():
    registry = load_registry()
    assert registry.families
    for family in registry.families.values():
        assert family.familywise_alpha == 0.01
        assert family.member_alpha * len(family.members) == pytest.approx(0.01)
        assert set(family.members) <= set(registry.by_name)
    assert registry.marker_alpha("CD11b") == pytest.approx(0.001)
    assert registry.marker_alpha("SST") == pytest.approx(0.01 / 3)
    # INS participates in two fixed families and takes the stricter allocation.
    assert registry.marker_alpha("INS") == pytest.approx(0.01 / 7)


def test_measurement_and_spillover_policy_are_explicit():
    registry = load_registry()
    for marker in registry.markers:
        assert marker.compartment in registry.allowed_compartments
        assert marker.spillover_policy in registry.allowed_spillover_policies
        assert marker.measurement_column == f"{marker.compartment}__{marker.name}"
        if marker.spillover_policy == "uncompensated":
            assert marker.compartment == "Nucleus"
    for marker in registry.active_markers:
        assert marker.compartment == "Cell"
        assert marker.spillover_policy == "redsea_corrected"
        assert marker.mutually_exclusive


def test_reference_graph_has_registered_nodes_and_allows_reciprocal_pairs():
    registry = load_registry()
    graph = registry.reference_graph
    assert graph["E_cadherin"] == "CD31"
    assert graph["CD31"] == "E_cadherin"
    for target, reference in graph.items():
        assert target != reference
        assert reference in registry.by_name
        assert target in registry.targets_for_reference(reference)


def test_registry_rejects_stale_cd11b_or_broad_sst():
    payload = json.loads(DEFAULT_REGISTRY_PATH.read_text())
    stale_cd11b = copy.deepcopy(payload)
    next(
        row for row in stale_cd11b["markers"] if row["name"] == "CD11b"
    )["reference"] = "EpCAM"
    with pytest.raises(ValueError, match="CD11b must use E_cadherin"):
        registry_from_dict(stale_cd11b)

    broad_sst = copy.deepcopy(payload)
    sst = next(row for row in broad_sst["markers"] if row["name"] == "SST")
    sst["uses"] = ["broad", "subtype"]
    sst["broad_compartment"] = "Epithelial"
    with pytest.raises(ValueError, match="SST must be active subtype-only"):
        registry_from_dict(broad_sst)


def test_registry_rejects_one_sided_family_membership():
    payload = json.loads(DEFAULT_REGISTRY_PATH.read_text())
    cd3e = next(row for row in payload["markers"] if row["name"] == "CD3e")
    cd3e["families"] = []
    with pytest.raises(ValueError, match="active calibration requires an error family"):
        registry_from_dict(payload)
