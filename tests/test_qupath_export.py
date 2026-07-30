import pandas as pd
import pytest

from phenocycler.qupath_export import export_assignment_frame


def test_export_preserves_uncertainty_and_compatibility_columns(tmp_path):
    assignments = pd.DataFrame(
        {
            "object_id": ["a", "b"],
            "donor_id": ["D", "D"],
            "broad_type": ["Immune", "Epithelial"],
            "specific_type": ["T", "Endocrine-delta"],
            "assignment_status": ["anchor", "ambiguous"],
            "confidence": [0.99, 0.55],
            "best_broad_type": ["Immune", "Epithelial"],
            "best_specific_type": ["T", "Endocrine-delta"],
            "broad_assignment_status": ["anchor", "anchor"],
            "specific_assignment_status": ["anchor", "ambiguous"],
            "ranked_broad_probabilities": [
                '[["Immune",0.99]]',
                '[["Epithelial",0.55]]',
            ],
            "ranked_specific_probabilities": [
                '[["T",0.99]]',
                '[["Endocrine-delta",0.55]]',
            ],
            "ranked_probabilities": ['[["Immune",0.99]]', '[["Epithelial",0.55]]'],
        }
    )
    context = pd.DataFrame({"object_id": ["a", "b"], "image": ["image.qptiff"] * 2})
    paths = export_assignment_frame(assignments, context, tmp_path)
    got = pd.read_csv(paths[0])
    assert list(got.assignment_status) == ["anchor", "ambiguous"]
    assert list(got.compartment) == list(got.broad_type)
    assert list(got.cell_type) == list(got.specific_type)


def test_export_rejects_missing_context(tmp_path):
    assignments = pd.DataFrame(
        {
            "object_id": ["a"],
            "donor_id": ["D"],
            "broad_type": ["Other"],
            "specific_type": ["Other"],
            "assignment_status": ["other"],
            "confidence": [1.0],
            "best_broad_type": ["Immune"],
            "best_specific_type": ["Other"],
            "broad_assignment_status": ["other"],
            "specific_assignment_status": ["other"],
            "ranked_broad_probabilities": ['[["Immune",0.4]]'],
            "ranked_specific_probabilities": ['[["Immune",0.4]]'],
            "ranked_probabilities": ['[["Immune",0.4]]'],
        }
    )
    with pytest.raises(ValueError, match="no image context"):
        export_assignment_frame(
            assignments,
            pd.DataFrame({"object_id": ["b"], "image": ["x"]}),
            tmp_path,
        )
