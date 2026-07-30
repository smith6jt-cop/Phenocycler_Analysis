"""Tests for the registry-driven simultaneous/hierarchical typing track."""

from __future__ import annotations

import copy
import json
import numpy as np
import pandas as pd
import pytest

from phenocycler import hierarchical_typing as ht


REGISTRY = {
    "marker_roles": {
        "CD3e": "identity",
        "CD68": "identity",
        "E_cadherin": "identity",
        "Pan_Cytokeratin": "identity",
        "Vimentin": "identity",
        "INS": "subtype_only",
        "SST": "identity",  # the module must override this unsafe declaration
        "Ki67": "process",
    },
    "broad": {
        "Immune": {
            "anchors": {"CD3e": 4.0},
            "support": {"CD68": 3.0},
            "exclusions": {"E_cadherin": 2.0},
            "required": ["CD3e", "CD68"],
        },
        "Epithelial": {
            "anchors": {"E_cadherin": 4.0},
            # SST and Ki67 are deliberately mislisted: role filtering must keep
            # both out of a broad identity model.
            "support": {"Pan_Cytokeratin": 2.0, "SST": 8.0, "Ki67": 8.0},
            "exclusions": {"CD3e": 2.0},
            "required": ["E_cadherin", "Pan_Cytokeratin", "SST", "Ki67"],
        },
        "Mesenchymal": {
            "anchors": {"Vimentin": 4.0},
            "exclusions": {"E_cadherin": 2.0},
            "required": ["Vimentin"],
        },
    },
    "subtypes": {
        "T_cell": {
            "parent": "Immune",
            "anchors": {"CD3e": 4.0},
            "required": ["CD3e"],
        },
        "Myeloid": {
            "parent": "Immune",
            "anchors": {"CD68": 4.0},
            "required": ["CD68"],
        },
        "Beta": {
            "parent": "Epithelial",
            "anchors": {"INS": 4.0},
            "required": ["INS"],
        },
        "Delta": {
            "parent": "Epithelial",
            "anchors": {"SST": 4.0},
            "required": ["SST"],
        },
    },
}


def _evidence(cells: dict[str, dict[str, object]], *, donor: str = "D1") -> pd.DataFrame:
    """Build long evidence.

    A marker value may be a probability, a state string, or ``(probability,
    state)``.
    """
    rows = []
    for object_id, markers in cells.items():
        for marker, value in markers.items():
            row = {"donor_id": donor, "object_id": object_id, "marker": marker}
            if isinstance(value, tuple):
                row["evidence_probability"], row["state"] = value
            elif isinstance(value, str):
                row["state"] = value
            else:
                row["evidence_probability"] = value
            rows.append(row)
    return pd.DataFrame(rows)


def _broad_negatives(**overrides):
    markers = {
        "CD3e": "negative",
        "CD68": "negative",
        "E_cadherin": "negative",
        "Pan_Cytokeratin": "negative",
        "Vimentin": "negative",
    }
    markers.update(overrides)
    return markers


def test_evidence_sources_normalize_to_one_probability_contract():
    raw = pd.DataFrame(
        [
            {"donor": "D", "object": "direct", "marker": "M", "evidence_probability": 0.8},
            {"donor": "D", "object": "tail", "marker": "M", "tail_probability": 0.01},
            {"donor": "D", "object": "ratio", "marker": "M", "log2_ratio": 2.0},
            {"donor": "D", "object": "state", "marker": "M", "state": "positive"},
            {
                "donor": "D",
                "object": "missing",
                "marker": "M",
                "evidence_probability": 0.99,
                "state": "unavailable",
            },
            {
                "donor": "D",
                "object": "uncertain",
                "marker": "M",
                "evidence_probability": 0.7,
                "state": "indeterminate",
            },
        ]
    )
    out = ht.prepare_cell_marker_evidence(raw).set_index("object_id")
    assert out.loc["direct", "evidence_probability"] == pytest.approx(0.8)
    assert out.loc["tail", "evidence_probability"] == pytest.approx(0.99)
    assert out.loc["ratio", "evidence_probability"] == pytest.approx(0.8)
    assert out.loc["state", "evidence_probability"] == 1.0
    assert np.isnan(out.loc["missing", "evidence_probability"])
    assert not out.loc["missing", "evidence_available"]
    assert out.loc["uncertain", "evidence_probability"] == pytest.approx(0.7)
    assert not out.loc["uncertain", "model_valid"]


def test_categorical_calibration_states_accept_missing_values():
    raw = pd.DataFrame(
        {
            "donor_id": ["D"] * 5,
            "object_id": ["p", "n", "i", "u", "missing"],
            "marker": ["M"] * 5,
            "state": pd.Categorical(
                [
                    "positive",
                    "negative",
                    "indeterminate",
                    "unavailable",
                    None,
                ],
                categories=[
                    "positive",
                    "negative",
                    "indeterminate",
                    "unavailable",
                ],
            ),
        }
    )

    out = ht.prepare_cell_marker_evidence(raw).set_index("object_id")

    assert out.loc["p", "state"] == "positive"
    assert out.loc["n", "state"] == "negative"
    assert out.loc["i", "state"] == "indeterminate"
    assert out.loc["u", "state"] == "unavailable"
    assert out.loc["missing", "state"] == ""
    assert np.isnan(out.loc["missing", "evidence_probability"])
    assert not out.loc["missing", "evidence_available"]


def test_wide_production_evidence_matches_long_assignment_without_expansion():
    long = _evidence(
        {
            "immune": _broad_negatives(CD3e=0.99),
            "other": _broad_negatives(),
            "missing": _broad_negatives(CD68="unavailable"),
        }
    )
    # Construct the compact production columns explicitly.
    production = long[["donor_id", "object_id"]].drop_duplicates().reset_index(drop=True)
    for marker in sorted(long["marker"].unique()):
        rows = long[long["marker"] == marker].set_index(["donor_id", "object_id"])
        key = pd.MultiIndex.from_frame(production[["donor_id", "object_id"]])
        if "evidence_probability" in rows:
            production[f"{marker}__expression_probability"] = rows[
                "evidence_probability"
            ].reindex(key).to_numpy()
        production[f"{marker}__state"] = rows["state"].reindex(key).fillna("").to_numpy()

    expected = ht.assign_broad_types(long, REGISTRY, stability=1.0).set_index("object_id")
    observed = ht.assign_broad_types(production, REGISTRY, stability=1.0).set_index("object_id")
    assert observed["broad_type"].to_dict() == expected["broad_type"].to_dict()
    assert (
        observed["assignment_status"].to_dict()
        == expected["assignment_status"].to_dict()
    )


def test_registry_enforces_identity_only_broad_and_sst_subtype_only():
    registry = ht.TypingRegistry.from_mapping(REGISTRY)
    assert registry.marker_role("SST") == ht.SUBTYPE_ONLY_ROLE
    broad = ht.classifier_from_registry(registry)
    assert "SST" not in broad.markers
    assert "Ki67" not in broad.markers
    epithelial = ht.classifier_from_registry(registry, level="subtype", parent="Epithelial")
    assert "SST" in epithelial.markers
    assert "INS" in epithelial.markers
    assert "Ki67" not in epithelial.markers


def test_calibration_registry_adapter_reuses_active_marker_membership():
    from phenocycler.marker_registry import load_registry

    source = load_registry()
    registry = ht.TypingRegistry.from_marker_registry(
        source,
        anchors={
            "Immune": ["CD3e"],
            "Endothelial": ["CD31"],
            "Epithelial": ["E_cadherin"],
            "Neural": ["B3TUBB"],
            "Mesenchymal": ["Vimentin"],
        },
    )
    expected = {}
    for marker in source.markers:
        if marker.kind == "type" and marker.active and "broad" in marker.uses:
            expected.setdefault(marker.broad_compartment, set()).add(marker.name)
    actual = {
        rule.name: set(rule.anchor_markers + rule.support_markers)
        for rule in registry.broad_rules
    }
    assert actual == expected
    assert "SST" not in set().union(*actual.values())
    assert registry.marker_role("SST") == ht.SUBTYPE_ONLY_ROLE
    assert registry.marker_role("Ki67") == ht.PROCESS_ROLE


def test_default_calibration_adapter_loads_versioned_biological_overlay():
    from phenocycler.marker_registry import load_registry

    source = load_registry()
    overlay = ht.load_typing_rules()
    registry = ht.TypingRegistry.from_marker_registry(source)

    assert registry.rules_version == "hierarchical-typing-rules-v1"
    assert registry.rules_fingerprint == overlay.fingerprint
    assert registry.marker_registry_fingerprint == source.fingerprint
    assert len(registry.rules_fingerprint) == 64

    immune = next(rule for rule in registry.broad_rules if rule.name == "Immune")
    assert "CD11b" in immune.support_markers
    assert source.marker("CD11b").reference == "E_cadherin"

    broad_positive = {
        marker
        for rule in registry.broad_rules
        for marker in rule.anchor_markers + rule.support_markers
    }
    assert "SST" not in broad_positive
    sst_subtypes = [
        rule.name
        for rule in registry.subtype_rules
        if "SST" in rule.anchor_markers + rule.support_markers + rule.exclusion_markers
    ]
    assert sst_subtypes == ["Delta"]

    state_markers = {marker.name for marker in source.state_markers}
    every_model_marker = {
        marker
        for rule in registry.broad_rules + registry.subtype_rules
        for marker in rule.model_markers
    }
    assert every_model_marker.isdisjoint(state_markers)


def test_typing_rules_fingerprint_is_canonical_and_validation_fails_loud():
    from phenocycler.marker_registry import load_registry

    source = load_registry()
    payload = json.loads(ht.DEFAULT_TYPING_RULES_PATH.read_text())
    canonical_copy = json.loads(json.dumps(payload, sort_keys=True))
    assert (
        ht.typing_rules_from_dict(payload).fingerprint
        == ht.typing_rules_from_dict(canonical_copy).fingerprint
    )

    wrong_version = copy.deepcopy(payload)
    wrong_version["marker_registry_version"] = "some-other-registry"
    with pytest.raises(ValueError, match="overlay expects marker registry"):
        ht.validate_typing_rules(ht.typing_rules_from_dict(wrong_version), source)

    process_leak = copy.deepcopy(payload)
    process_leak["broad"]["Immune"]["support"].append("Ki67")
    with pytest.raises(ValueError, match="process/state marker 'Ki67'"):
        ht.validate_typing_rules(ht.typing_rules_from_dict(process_leak), source)

    sst_leak = copy.deepcopy(payload)
    sst_leak["subtypes"]["Beta"]["support"].append("SST")
    with pytest.raises(ValueError, match="SST must appear exactly once"):
        ht.validate_typing_rules(ht.typing_rules_from_dict(sst_leak), source)


def test_explicit_typing_rules_cannot_mix_with_ad_hoc_overrides():
    from phenocycler.marker_registry import load_registry

    with pytest.raises(ValueError, match="cannot be combined"):
        ht.TypingRegistry.from_marker_registry(
            load_registry(),
            typing_rules=ht.load_typing_rules(),
            anchors={"Immune": ["CD3e"]},
        )


def test_balanced_seed_weights_equalize_class_then_donor():
    seeds = pd.DataFrame(
        {
            "donor_id": ["A", "A", "B", "C", "C", "C", "C"],
            "label": ["Immune", "Immune", "Immune", "Epithelial", "Epithelial", "Epithelial", "Epithelial"],
        }
    )
    weights = ht.balanced_seed_weights(seeds)
    audit = seeds.assign(weight=weights)
    by_class = audit.groupby("label").weight.sum()
    assert by_class["Immune"] == pytest.approx(0.5)
    assert by_class["Epithelial"] == pytest.approx(0.5)
    by_donor = audit[audit.label == "Immune"].groupby("donor_id").weight.sum()
    assert by_donor["A"] == pytest.approx(0.25)
    assert by_donor["B"] == pytest.approx(0.25)
    assert weights.sum() == pytest.approx(1.0)


def _training_fixture():
    cells = {}
    seed_rows = []
    for donor in ("D1", "D2", "D3"):
        immune = f"{donor}_immune"
        epithelial = f"{donor}_epi"
        mesenchymal = f"{donor}_mes"
        cells[(donor, immune)] = _broad_negatives(CD3e=0.99, CD68=0.9)
        cells[(donor, epithelial)] = _broad_negatives(
            E_cadherin=0.99, Pan_Cytokeratin=0.9
        )
        cells[(donor, mesenchymal)] = _broad_negatives(Vimentin=0.99)
        seed_rows.extend(
            [
                {"donor_id": donor, "object_id": immune, "label": "Immune"},
                {"donor_id": donor, "object_id": epithelial, "label": "Epithelial"},
                {"donor_id": donor, "object_id": mesenchymal, "label": "Mesenchymal"},
            ]
        )
    rows = []
    for (donor, object_id), markers in cells.items():
        part = _evidence({object_id: markers}, donor=donor)
        rows.append(part)
    return pd.concat(rows, ignore_index=True), pd.DataFrame(seed_rows)


def test_fitted_multinomial_is_additive_sign_constrained_and_uniform_prior():
    evidence, seeds = _training_fixture()
    model = ht.fit_constrained_classifier(evidence, REGISTRY, seeds, l2=0.01)
    coefficients = model.coefficient_table()
    assert model.class_priors == {
        "Immune": pytest.approx(1 / 3),
        "Epithelial": pytest.approx(1 / 3),
        "Mesenchymal": pytest.approx(1 / 3),
    }
    assert coefficients.loc["Immune", "CD3e"] >= 0
    assert coefficients.loc["Immune", "CD68"] >= 0
    assert coefficients.loc["Immune", "E_cadherin"] <= 0
    assert coefficients.loc["Epithelial", "E_cadherin"] >= 0
    assert coefficients.loc["Epithelial", "CD3e"] <= 0
    # A marker not curated for a class is fixed exactly at zero.
    assert coefficients.loc["Mesenchymal", "CD68"] == 0.0
    # With no marker evidence and no intercept, priors are exactly uniform.
    blank = _evidence({"blank": {"Ki67": 0.99}})
    probability = model.predict_proba(blank)
    assert probability.loc[0, list(model.classes)].tolist() == pytest.approx([1 / 3] * 3)


def test_broad_assignment_contract_covers_every_authoritative_state():
    evidence = _evidence(
        {
            "anchor": _broad_negatives(CD3e=0.99),
            "anchor_conflict": _broad_negatives(CD3e=0.99, E_cadherin=0.99),
            "inferred": _broad_negatives(CD68=0.99),
            "ambiguous": _broad_negatives(CD68=0.35),
            "other": _broad_negatives(),
            "unavailable": _broad_negatives(Pan_Cytokeratin="unavailable"),
            "process_only": {"Ki67": 0.99},
        }
    )
    out = ht.assign_broad_types(evidence, REGISTRY, stability=0.95).set_index("object_id")

    assert out.loc["anchor", "broad_label"] == "Immune"
    assert out.loc["anchor", "broad_assignment_status"] == ht.AUTHORITATIVE_ANCHOR
    assert out.loc["anchor_conflict", "broad_label"] == ht.AMBIGUOUS_LABEL
    assert out.loc["anchor_conflict", "broad_assignment_status"] == ht.AMBIGUOUS
    assert out.loc["inferred", "broad_label"] == "Immune"
    assert out.loc["inferred", "broad_assignment_status"] == ht.INFERRED
    assert out.loc["ambiguous", "broad_assignment_status"] == ht.AMBIGUOUS
    assert out.loc["other", "broad_label"] == ht.OTHER_LABEL
    assert out.loc["other", "broad_assignment_status"] == ht.OTHER
    assert out.loc["unavailable", "broad_label"] == ht.UNAVAILABLE_LABEL
    assert out.loc["unavailable", "broad_assignment_status"] == ht.UNAVAILABLE
    assert out.loc["process_only", "broad_assignment_status"] == ht.UNAVAILABLE

    # A candidate is an audit field, not an authoritative assignment. It is
    # retained for ambiguous/other/incomplete cells whenever identity evidence
    # exists, and absent for a process-only cell.
    for object_id in ("anchor", "anchor_conflict", "inferred", "ambiguous", "other", "unavailable"):
        assert out.loc[object_id, "broad_best_candidate"] in {
            "Immune",
            "Epithelial",
            "Mesenchymal",
        }
    assert out.loc["process_only", "broad_best_candidate"] == ""


def test_negative_calibration_state_cannot_become_anchor_from_rank_probability():
    evidence = _evidence(
        {
            "subthreshold": {
                **_broad_negatives(),
                "CD3e": (0.999, "negative"),
            }
        }
    )
    out = ht.assign_broad_types(evidence, REGISTRY, stability=1.0).iloc[0]
    assert out["broad_assignment_status"] == ht.OTHER
    assert out["broad_label"] == ht.OTHER_LABEL
    assert out["broad_authoritative_markers"] == ""


def test_inference_uses_exact_probability_margin_and_stability_thresholds():
    evidence = _evidence(
        {
            "unstable": _broad_negatives(CD68=0.99),
            "stable": _broad_negatives(CD68=0.99),
        }
    )
    stability = pd.DataFrame(
        {
            "donor_id": ["D1", "D1"],
            "object_id": ["unstable", "stable"],
            "stability": [0.899, 0.90],
        }
    )
    out = ht.assign_broad_types(evidence, REGISTRY, stability=stability).set_index("object_id")
    assert out.loc["unstable", "broad_assignment_status"] == ht.AMBIGUOUS
    assert out.loc["stable", "broad_assignment_status"] == ht.INFERRED
    assert out.loc["stable", "broad_best_probability"] >= 0.80
    assert out.loc["stable", "broad_margin"] >= 0.25


def test_other_requires_every_required_model_valid_and_negative():
    evidence = _evidence(
        {
            "complete": _broad_negatives(),
            "missing_one": _broad_negatives(CD68="unavailable"),
            "borderline": _broad_negatives(CD68=0.3),
        }
    )
    out = ht.assign_broad_types(evidence, REGISTRY, stability=1.0).set_index("object_id")
    assert out.loc["complete", "broad_assignment_status"] == ht.OTHER
    assert out.loc["complete", "broad_all_required_valid"]
    assert out.loc["complete", "broad_all_required_negative"]
    assert out.loc["missing_one", "broad_assignment_status"] == ht.UNAVAILABLE
    assert not out.loc["missing_one", "broad_all_required_valid"]
    assert out.loc["borderline", "broad_assignment_status"] == ht.AMBIGUOUS
    assert out.loc["borderline", "broad_all_required_valid"]
    assert not out.loc["borderline", "broad_all_required_negative"]


def test_sst_is_subtype_only_and_subtypes_cannot_escape_parent():
    evidence = _evidence(
        {
            "delta": {
                **_broad_negatives(E_cadherin=0.99),
                "INS": "negative",
                "SST": 0.99,
            },
            "immune_with_sst": {
                **_broad_negatives(CD3e=0.99),
                "INS": "negative",
                "SST": 0.99,
            },
            "epi_unclassified": {
                **_broad_negatives(E_cadherin=0.99),
                "INS": "negative",
                "SST": "negative",
            },
            "terminal_parent": {
                **_broad_negatives(Vimentin=0.99),
                "SST": 0.99,
            },
            "sst_only": {"SST": 0.99},
        }
    )
    out = ht.assign_hierarchical_types(
        evidence,
        REGISTRY,
        broad_stability=1.0,
        subtype_stability={"Immune": 1.0, "Epithelial": 1.0},
    ).set_index("object_id")

    assert out.loc["delta", "broad_label"] == "Epithelial"
    assert out.loc["delta", "cell_type"] == "Delta"
    assert out.loc["delta", "subtype_assignment_status"] == ht.AUTHORITATIVE_ANCHOR

    assert out.loc["immune_with_sst", "broad_label"] == "Immune"
    assert out.loc["immune_with_sst", "cell_type"] == "T_cell"
    assert out.loc["immune_with_sst", "cell_type"] != "Delta"

    assert out.loc["epi_unclassified", "broad_label"] == "Epithelial"
    assert out.loc["epi_unclassified", "cell_type"] == "Epithelial-unclassified"
    assert out.loc["epi_unclassified", "subtype_assignment_status"] == ht.OTHER

    assert out.loc["terminal_parent", "cell_type"] == "Mesenchymal"
    assert out.loc["terminal_parent", "subtype_assignment_status"] == "not_applicable"

    assert out.loc["sst_only", "broad_label"] == ht.UNAVAILABLE_LABEL
    assert out.loc["sst_only", "broad_best_candidate"] == ""
    assert out.loc["sst_only", "cell_type"] == ht.UNAVAILABLE_LABEL

    for column in (
        "broad_type",
        "specific_type",
        "assignment_status",
        "confidence",
        "best_broad_type",
        "best_specific_type",
        "specific_assignment_status",
        "ranked_probabilities",
        "ranked_broad_probabilities",
        "ranked_specific_probabilities",
    ):
        assert column in out.columns
    assert out.loc["delta", "broad_type"] == "Epithelial"
    assert out.loc["delta", "specific_type"] == "Delta"
    assert out.loc["delta", "best_specific_type"] == "Delta"
    assert out.loc["delta", "assignment_status"] == "anchor"
    ranking = json.loads(out.loc["delta", "ranked_probabilities"])
    assert ranking[0][0] == out.loc["delta", "best_specific_type"]
    assert sum(probability for _, probability in ranking) == pytest.approx(1.0)


def test_grouped_donor_bootstrap_returns_cell_level_candidate_stability():
    evidence, seeds = _training_fixture()
    result = ht.grouped_donor_bootstrap_stability(
        evidence,
        REGISTRY,
        seeds,
        n_bootstrap=5,
        random_state=7,
        maxiter=100,
    )
    assert len(result) == evidence[["donor_id", "object_id"]].drop_duplicates().shape[0]
    assert (result["bootstrap_replicates"] == 5).all()
    assert result["stability"].between(0, 1).all()
    assert result["best_candidate"].isin({"Immune", "Epithelial", "Mesenchymal"}).all()
    assert not result.duplicated(["donor_id", "object_id"]).any()
