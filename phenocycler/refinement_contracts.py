"""Immutable contracts for donor-disjoint probabilistic typing development.

These contracts deliberately live outside the production assignment algorithm.
They describe which donors and independently adjudicated labels may be used to
fit, calibrate, validate, or lock a candidate model.  Existing QuPath labels and
the historical ``gold_standard`` bundles are audit-only and cannot satisfy the
development-label purpose asserted here.
"""

from __future__ import annotations

import hashlib
import hmac
import json
import os
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np
import pandas as pd

from .artifacts import (
    DatasetSnapshot,
    RunManifest,
    StageManifest,
    dataframe_content_sha256,
    sha256_json,
    validate_dataset_snapshot_current,
)

SPLIT_SCHEMA_VERSION = 2
LABEL_SCHEMA_VERSION = 2
REVIEW_LABEL_SCHEMA_VERSION = 5
EVIDENCE_SCHEMA_VERSION = 1
DEVELOPMENT_LABEL_PURPOSE = "probabilistic_typing_development"
DEVELOPMENT_LABEL_SOURCE = "independent_blinded_review"
DEVELOPMENT_EVIDENCE_PURPOSE = "probabilistic_typing_training_evidence"
REVIEW_LABEL_PURPOSES = {
    "development": "probabilistic_typing_development_review",
    "calibration": "probabilistic_typing_calibration_review",
    "validation": "probabilistic_typing_validation_review",
    "locked_test": "probabilistic_typing_locked_test_review",
}
SPLIT_NAMES = (
    "development",
    "calibration",
    "validation",
    "locked_test",
)
REVIEW_SPLIT_NAMES = ("calibration", "validation", "locked_test")
REVIEW_SAMPLING_METHOD = "hmac_sha256_candidate_bound_v2"
DEFAULT_REVIEW_N_PER_STRATUM = 40
LABEL_COLUMNS = (
    "donor_id",
    "object_id",
    "level",
    "parent",
    "label",
)
LABEL_GOVERNANCE_COLUMNS = (
    "evidence_sufficient",
    "audit_only",
    "label_purpose",
    "label_source",
    "reviewer_1",
    "reviewer_2",
    "adjudicated",
    "root_cause",
)


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> Path:
    """Write a new JSON artifact atomically and never overwrite an artifact."""

    path = path.expanduser().resolve()
    if path.exists():
        raise FileExistsError(f"immutable artifact already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        # Linking a fully written same-filesystem temporary into its final name
        # is an atomic create-if-absent operation. Unlike check-then-replace,
        # concurrent writers cannot overwrite an immutable artifact.
        os.link(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return path


def _strict_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"JSON contains duplicate key {key!r}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"JSON contains non-finite number {value!r}")


def _read_strict_json(path: Path | str) -> Mapping[str, Any]:
    with Path(path).expanduser().open(encoding="utf-8") as handle:
        value = json.load(
            handle,
            object_pairs_hook=_strict_json_object,
            parse_constant=_reject_json_constant,
        )
    if not isinstance(value, Mapping):
        raise TypeError("JSON artifact must contain one object")
    return value


def _normalise_donors(values: Sequence[object], *, split: str) -> tuple[str, ...]:
    if isinstance(values, (str, bytes)):
        raise TypeError(f"{split} donors must be a sequence, not a string")
    donors = tuple(str(value).strip() for value in values)
    if any(not donor for donor in donors):
        raise ValueError(f"{split} donors must be non-empty strings")
    if len(set(donors)) != len(donors):
        raise ValueError(f"{split} contains duplicate donors")
    return tuple(sorted(donors))


@dataclass(frozen=True)
class DonorSplitManifest:
    """A content-addressed, mutually exclusive donor split.

    Development donors fit coefficients and bootstrap ensembles. Calibration
    donors may fit probability calibration and acceptance thresholds. Validation
    and locked-test donors are recorded for later evaluation and are never read
    by the fitting API.
    """

    split_version: str
    source_run_id: str
    cohort_content_id: str
    development: tuple[str, ...]
    calibration: tuple[str, ...]
    validation: tuple[str, ...] = ()
    locked_test: tuple[str, ...] = ()
    review_n_per_stratum: Mapping[str, int] = field(
        default_factory=lambda: {
            split: DEFAULT_REVIEW_N_PER_STRATUM for split in REVIEW_SPLIT_NAMES
        }
    )
    review_blinding_key_sha256: str = ""
    threshold_selection_policy_content_id: str = ""
    threshold_grid_content_id: str = ""
    promotion_policy_content_id: str = ""
    release_signing_key_sha256: str = ""
    schema_version: int = SPLIT_SCHEMA_VERSION
    content_id: str = ""

    def __post_init__(self) -> None:
        if int(self.schema_version) != SPLIT_SCHEMA_VERSION:
            raise ValueError(f"unsupported donor-split schema_version {self.schema_version!r}")
        for field_name in ("split_version", "source_run_id", "cohort_content_id"):
            if not str(getattr(self, field_name)).strip():
                raise ValueError(f"{field_name} must be a non-empty string")
        for name in SPLIT_NAMES:
            normalised = _normalise_donors(getattr(self, name), split=name)
            object.__setattr__(self, name, normalised)
        if not isinstance(self.review_n_per_stratum, Mapping):
            raise TypeError("review_n_per_stratum must be a mapping")
        unknown_review_splits = sorted(
            set(self.review_n_per_stratum) - set(REVIEW_SPLIT_NAMES)
        )
        missing_review_splits = sorted(
            set(REVIEW_SPLIT_NAMES) - set(self.review_n_per_stratum)
        )
        if unknown_review_splits or missing_review_splits:
            raise ValueError(
                "review_n_per_stratum must cover exactly calibration, validation, "
                f"and locked_test; missing={missing_review_splits}, "
                f"unknown={unknown_review_splits}"
            )
        review_n: dict[str, int] = {}
        for split in REVIEW_SPLIT_NAMES:
            value = self.review_n_per_stratum[split]
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise TypeError(f"review_n_per_stratum[{split!r}] must be an integer")
            if int(value) < 1:
                raise ValueError(f"review_n_per_stratum[{split!r}] must be positive")
            review_n[split] = int(value)
        object.__setattr__(
            self, "review_n_per_stratum", MappingProxyType(review_n)
        )
        commitment = str(self.review_blinding_key_sha256).strip().lower()
        if commitment and (
            len(commitment) != 64
            or any(character not in "0123456789abcdef" for character in commitment)
        ):
            raise ValueError(
                "review_blinding_key_sha256 must be blank or a lowercase SHA-256 digest"
            )
        object.__setattr__(self, "review_blinding_key_sha256", commitment)
        for field_name in (
            "threshold_selection_policy_content_id",
            "threshold_grid_content_id",
            "promotion_policy_content_id",
            "release_signing_key_sha256",
        ):
            digest = str(getattr(self, field_name)).strip().lower()
            if digest and (
                len(digest) != 64
                or any(character not in "0123456789abcdef" for character in digest)
            ):
                raise ValueError(f"{field_name} must be blank or a lowercase SHA-256 digest")
            object.__setattr__(self, field_name, digest)
        if not self.development:
            raise ValueError("development must contain at least one donor")
        if not self.calibration:
            raise ValueError("calibration must contain at least one donor")
        owner: dict[str, str] = {}
        overlap: dict[str, list[str]] = {}
        for split in SPLIT_NAMES:
            for donor in getattr(self, split):
                if donor in owner:
                    overlap.setdefault(donor, [owner[donor]]).append(split)
                else:
                    owner[donor] = split
        if overlap:
            raise ValueError(f"donor splits overlap: {overlap}")
        expected = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected:
            raise ValueError("donor split content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    @property
    def donors(self) -> tuple[str, ...]:
        return tuple(donor for split in SPLIT_NAMES for donor in getattr(self, split))

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": int(self.schema_version),
            "split_version": self.split_version,
            "source_run_id": self.source_run_id,
            "cohort_content_id": self.cohort_content_id,
            "splits": {name: list(getattr(self, name)) for name in SPLIT_NAMES},
            "review_sampling": {
                "method": REVIEW_SAMPLING_METHOD,
                "n_per_stratum": {
                    name: self.review_n_per_stratum[name]
                    for name in REVIEW_SPLIT_NAMES
                },
                "blinding_key_sha256": self.review_blinding_key_sha256,
            },
            "release_governance": {
                "threshold_selection_policy_content_id": (
                    self.threshold_selection_policy_content_id
                ),
                "threshold_grid_content_id": self.threshold_grid_content_id,
                "promotion_policy_content_id": self.promotion_policy_content_id,
                "signing_key_sha256": self.release_signing_key_sha256,
            },
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    def split_for(self, donor_id: object) -> str:
        donor = str(donor_id)
        for split in SPLIT_NAMES:
            if donor in getattr(self, split):
                return split
        raise KeyError(f"donor {donor!r} is absent from split manifest")

    def validate_cohort(self, expected_donors: Sequence[object]) -> None:
        expected = {str(value) for value in expected_donors}
        observed = set(self.donors)
        missing = sorted(expected - observed)
        unknown = sorted(observed - expected)
        if missing or unknown:
            raise ValueError(
                "donor split is not an exact cohort partition: "
                f"missing={missing}, unknown={unknown}"
            )

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> DonorSplitManifest:
        if not isinstance(value, Mapping):
            raise TypeError("donor split manifest must be a mapping")
        allowed = {
            "schema_version",
            "split_version",
            "source_run_id",
            "cohort_content_id",
            "splits",
            "review_sampling",
            "release_governance",
            "content_id",
        }
        unknown_fields = sorted(set(value) - allowed)
        if unknown_fields:
            raise ValueError(f"donor split manifest has unknown fields: {unknown_fields}")
        if not str(value.get("content_id", "")).strip():
            raise ValueError("donor split manifest is missing content_id")
        splits = value.get("splits")
        if not isinstance(splits, Mapping):
            raise ValueError("donor split manifest is missing a splits mapping")
        unknown = sorted(set(splits) - set(SPLIT_NAMES))
        if unknown:
            raise ValueError(f"unknown donor split names: {unknown}")
        review_sampling = value.get("review_sampling")
        if not isinstance(review_sampling, Mapping):
            raise ValueError("donor split manifest is missing review_sampling")
        review_fields = {"method", "n_per_stratum", "blinding_key_sha256"}
        if set(review_sampling) != review_fields:
            raise ValueError(
                "review_sampling must contain exactly method, n_per_stratum, "
                "and blinding_key_sha256"
            )
        if review_sampling["method"] != REVIEW_SAMPLING_METHOD:
            raise ValueError("donor split uses an unsupported review sampling method")
        review_n = review_sampling["n_per_stratum"]
        if not isinstance(review_n, Mapping):
            raise TypeError("review_sampling.n_per_stratum must be a mapping")
        release_governance = value.get("release_governance")
        if not isinstance(release_governance, Mapping) or set(release_governance) != {
            "threshold_selection_policy_content_id",
            "threshold_grid_content_id",
            "promotion_policy_content_id",
            "signing_key_sha256",
        }:
            raise ValueError(
                "release_governance must contain exactly "
                "threshold-selection policy/grid IDs, promotion policy ID, "
                "and signing_key_sha256"
            )
        return cls(
            schema_version=int(value.get("schema_version", -1)),
            split_version=str(value.get("split_version", "")),
            source_run_id=str(value.get("source_run_id", "")),
            cohort_content_id=str(value.get("cohort_content_id", "")),
            development=tuple(splits.get("development", ())),
            calibration=tuple(splits.get("calibration", ())),
            validation=tuple(splits.get("validation", ())),
            locked_test=tuple(splits.get("locked_test", ())),
            review_n_per_stratum=review_n,
            review_blinding_key_sha256=review_sampling[
                "blinding_key_sha256"
            ],
            promotion_policy_content_id=release_governance[
                "promotion_policy_content_id"
            ],
            threshold_selection_policy_content_id=release_governance[
                "threshold_selection_policy_content_id"
            ],
            threshold_grid_content_id=release_governance[
                "threshold_grid_content_id"
            ],
            release_signing_key_sha256=release_governance[
                "signing_key_sha256"
            ],
            content_id=str(value.get("content_id", "")),
        )

    @classmethod
    def read_json(cls, path: Path | str) -> DonorSplitManifest:
        return cls.from_dict(_read_strict_json(path))


def _canonical_label_frame(labels: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(labels, pd.DataFrame):
        raise TypeError("labels must be a pandas DataFrame")
    work = labels.copy()
    if "parent" not in work:
        work["parent"] = ""
    missing = sorted(set((*LABEL_COLUMNS, *LABEL_GOVERNANCE_COLUMNS)) - set(work.columns))
    if missing:
        raise KeyError(f"development labels are missing columns: {missing}")
    for column in LABEL_COLUMNS:
        work[column] = work[column].astype("string").fillna("").str.strip()
    work["level"] = work["level"].str.lower()
    if work.empty:
        raise ValueError("development label bundle is empty")
    if (work[["donor_id", "object_id", "label"]] == "").any().any():
        raise ValueError("donor_id, object_id, and label must be non-empty")
    unknown_levels = sorted(set(work["level"]) - {"broad", "subtype"})
    if unknown_levels:
        raise ValueError(f"unknown label levels: {unknown_levels}")
    bad_parent = ((work["level"] == "broad") & (work["parent"] != "")) | (
        (work["level"] == "subtype") & (work["parent"] == "")
    )
    if bad_parent.any():
        raise ValueError("broad labels require blank parent; subtype labels require parent")
    key = ["donor_id", "object_id", "level", "parent"]
    if work.duplicated(key).any():
        examples = work.loc[work.duplicated(key, keep=False), key].head(5)
        raise ValueError(
            "development labels contain duplicate cell/level rows: "
            f"{examples.to_dict(orient='records')}"
        )
    subtype = work.loc[
        work["level"].eq("subtype"),
        ["donor_id", "object_id", "parent"],
    ]
    if not subtype.empty:
        broad = work.loc[
            work["level"].eq("broad"),
            ["donor_id", "object_id", "label"],
        ].rename(columns={"label": "broad_parent_label"})
        aligned = subtype.merge(
            broad,
            on=["donor_id", "object_id"],
            how="left",
            validate="m:1",
        )
        inconsistent = aligned["broad_parent_label"].isna() | aligned["broad_parent_label"].ne(
            aligned["parent"]
        )
        if inconsistent.any():
            examples = aligned.loc[
                inconsistent,
                ["donor_id", "object_id", "parent", "broad_parent_label"],
            ].head(5)
            raise ValueError(
                "subtype labels require a matching broad parent label: "
                f"{examples.to_dict(orient='records')}"
            )
    sufficient = work["evidence_sufficient"]
    if sufficient.isna().any() or not sufficient.isin([True, False, 0, 1]).all():
        raise ValueError("evidence_sufficient must contain explicit booleans")
    if not sufficient.astype(bool).all():
        raise ValueError("insufficient-evidence review rows cannot train typing")
    audit_only = work["audit_only"]
    if audit_only.isna().any() or not audit_only.isin([True, False, 0, 1]).all():
        raise ValueError("audit_only must contain explicit booleans")
    if audit_only.astype(bool).any():
        raise ValueError("audit-only labels cannot train probabilistic typing")
    purpose = work["label_purpose"].astype("string").fillna("").str.strip()
    if not purpose.eq(DEVELOPMENT_LABEL_PURPOSE).all():
        raise ValueError(
            f"every training row must declare label_purpose={DEVELOPMENT_LABEL_PURPOSE!r}"
        )
    label_source = work["label_source"].astype("string").fillna("").str.strip()
    if not label_source.eq(DEVELOPMENT_LABEL_SOURCE).all():
        raise ValueError(
            f"every training row must declare label_source={DEVELOPMENT_LABEL_SOURCE!r}"
        )
    reviewer_1 = work["reviewer_1"].astype("string").fillna("").str.strip()
    reviewer_2 = work["reviewer_2"].astype("string").fillna("").str.strip()
    root_cause = work["root_cause"].astype("string").fillna("").str.strip()
    if reviewer_1.eq("").any() or reviewer_2.eq("").any():
        raise ValueError("every training row requires two non-empty reviewers")
    if reviewer_1.eq(reviewer_2).any():
        raise ValueError("every training row requires two distinct reviewers")
    if root_cause.eq("").any():
        raise ValueError("every training row requires a non-empty root_cause")
    adjudicated = work["adjudicated"]
    if adjudicated.isna().any() or not adjudicated.isin([True, False, 0, 1]).all():
        raise ValueError("adjudicated must contain explicit booleans")
    if not adjudicated.astype(bool).all():
        raise ValueError("every training row must be adjudicated")
    work["evidence_sufficient"] = sufficient.astype(bool)
    work["audit_only"] = audit_only.astype(bool)
    work["label_purpose"] = purpose
    work["label_source"] = label_source
    work["reviewer_1"] = reviewer_1
    work["reviewer_2"] = reviewer_2
    work["adjudicated"] = adjudicated.astype(bool)
    work["root_cause"] = root_cause
    return work


def label_frame_fingerprint(labels: pd.DataFrame) -> str:
    """Content fingerprint independent of row order and unrelated columns."""

    work = _canonical_label_frame(labels)
    records = (
        work.loc[:, [*LABEL_COLUMNS, *LABEL_GOVERNANCE_COLUMNS]]
        .sort_values(list(LABEL_COLUMNS), kind="mergesort")
        .to_dict(orient="records")
    )
    return sha256_json({"schema_version": LABEL_SCHEMA_VERSION, "labels": records})


@dataclass(frozen=True)
class DevelopmentLabelProvenance:
    """Identity of a new label set explicitly authorized for model development."""

    bundle_id: str
    source_run_id: str
    split_manifest_content_id: str
    fingerprint: str
    purpose: str = DEVELOPMENT_LABEL_PURPOSE
    schema_version: int = LABEL_SCHEMA_VERSION
    content_id: str = ""

    def __post_init__(self) -> None:
        if int(self.schema_version) != LABEL_SCHEMA_VERSION:
            raise ValueError(
                f"unsupported development-label schema_version {self.schema_version!r}"
            )
        if self.purpose != DEVELOPMENT_LABEL_PURPOSE:
            raise ValueError(
                "labels are not explicitly authorized as a separate "
                f"{DEVELOPMENT_LABEL_PURPOSE!r} bundle"
            )
        for field_name in (
            "bundle_id",
            "source_run_id",
            "split_manifest_content_id",
            "fingerprint",
        ):
            if not str(getattr(self, field_name)).strip():
                raise ValueError(f"{field_name} must be a non-empty string")
        if len(self.split_manifest_content_id) != 64 or len(self.fingerprint) != 64:
            raise ValueError("development-label split/fingerprint must be SHA-256")
        expected = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected:
            raise ValueError("development-label provenance content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": int(self.schema_version),
            "bundle_id": self.bundle_id,
            "source_run_id": self.source_run_id,
            "split_manifest_content_id": self.split_manifest_content_id,
            "purpose": self.purpose,
            "fingerprint": self.fingerprint,
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_labels(
        cls,
        labels: pd.DataFrame,
        *,
        bundle_id: str,
        splits: DonorSplitManifest,
        purpose: str,
    ) -> DevelopmentLabelProvenance:
        return cls(
            bundle_id=str(bundle_id),
            source_run_id=splits.source_run_id,
            split_manifest_content_id=splits.content_id,
            purpose=str(purpose),
            fingerprint=label_frame_fingerprint(labels),
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> DevelopmentLabelProvenance:
        if not isinstance(value, Mapping):
            raise TypeError("development-label provenance must be a mapping")
        allowed = {
            "schema_version",
            "bundle_id",
            "source_run_id",
            "split_manifest_content_id",
            "purpose",
            "fingerprint",
            "content_id",
        }
        unknown = sorted(set(value) - allowed)
        if unknown:
            raise ValueError(f"development-label provenance has unknown fields: {unknown}")
        if not str(value.get("content_id", "")).strip():
            raise ValueError("development-label provenance is missing content_id")
        return cls(
            schema_version=int(value.get("schema_version", -1)),
            bundle_id=str(value.get("bundle_id", "")),
            source_run_id=str(value.get("source_run_id", "")),
            split_manifest_content_id=str(value.get("split_manifest_content_id", "")),
            purpose=str(value.get("purpose", "")),
            fingerprint=str(value.get("fingerprint", "")),
            content_id=str(value.get("content_id", "")),
        )

    @classmethod
    def read_json(cls, path: Path | str) -> DevelopmentLabelProvenance:
        return cls.from_dict(_read_strict_json(path))


def evidence_frame_fingerprint(evidence: pd.DataFrame) -> str:
    """Hash the exact labeled evidence used by fitting, independent of row order."""

    if not isinstance(evidence, pd.DataFrame):
        raise TypeError("evidence must be a pandas DataFrame")
    if evidence.columns.has_duplicates or any(
        not isinstance(column, str) or not column for column in evidence.columns
    ):
        raise ValueError("evidence columns must be unique non-empty strings")
    missing = sorted({"donor_id", "object_id"} - set(evidence.columns))
    if missing:
        raise KeyError(f"training evidence is missing keys: {missing}")
    work = evidence.copy()
    for column in ("donor_id", "object_id"):
        if work[column].isna().any():
            raise ValueError(f"training evidence {column} contains missing values")
        work[column] = work[column].astype(str)
        if work[column].eq("").any():
            raise ValueError(f"training evidence {column} contains blank values")
    key = ["donor_id", "object_id"]
    if "marker" in work:
        work["marker"] = work["marker"].astype("string").fillna("").astype(str)
        if work["marker"].eq("").any():
            raise ValueError("training evidence marker contains blank values")
        key.append("marker")
    if work.duplicated(key).any():
        raise ValueError(f"training evidence contains duplicate keys: {key}")
    ordered_columns = [*key, *sorted(set(work.columns) - set(key))]
    return dataframe_content_sha256(
        work.loc[:, ordered_columns],
        sort_by=key,
    )


@dataclass(frozen=True)
class DevelopmentEvidenceProvenance:
    """Identity of the immutable evidence artifact and fitted labeled projection."""

    source_run_id: str
    split_manifest_content_id: str
    source_run_manifest_content_id: str
    source_stage_manifest_content_id: str
    source_artifact_content_id: str
    labeled_frame_fingerprint: str
    feature_schema_version: str
    row_count: int
    purpose: str = DEVELOPMENT_EVIDENCE_PURPOSE
    schema_version: int = EVIDENCE_SCHEMA_VERSION
    content_id: str = ""

    def __post_init__(self) -> None:
        if int(self.schema_version) != EVIDENCE_SCHEMA_VERSION:
            raise ValueError(
                f"unsupported evidence-provenance schema_version {self.schema_version!r}"
            )
        if self.purpose != DEVELOPMENT_EVIDENCE_PURPOSE:
            raise ValueError(
                "training evidence is not explicitly authorized for probabilistic typing"
            )
        for field_name in (
            "source_run_id",
            "split_manifest_content_id",
            "source_run_manifest_content_id",
            "source_stage_manifest_content_id",
            "source_artifact_content_id",
            "labeled_frame_fingerprint",
            "feature_schema_version",
        ):
            if not str(getattr(self, field_name)).strip():
                raise ValueError(f"{field_name} must be a non-empty string")
        for field_name in (
            "split_manifest_content_id",
            "source_run_manifest_content_id",
            "source_stage_manifest_content_id",
            "source_artifact_content_id",
            "labeled_frame_fingerprint",
        ):
            if len(str(getattr(self, field_name))) != 64:
                raise ValueError(f"{field_name} must be SHA-256")
        if isinstance(self.row_count, bool) or int(self.row_count) < 1:
            raise ValueError("training evidence row_count must be a positive integer")
        object.__setattr__(self, "row_count", int(self.row_count))
        expected = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected:
            raise ValueError("evidence provenance content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": int(self.schema_version),
            "purpose": self.purpose,
            "source_run_id": self.source_run_id,
            "split_manifest_content_id": self.split_manifest_content_id,
            "source_run_manifest_content_id": self.source_run_manifest_content_id,
            "source_stage_manifest_content_id": self.source_stage_manifest_content_id,
            "source_artifact_content_id": self.source_artifact_content_id,
            "labeled_frame_fingerprint": self.labeled_frame_fingerprint,
            "feature_schema_version": self.feature_schema_version,
            "row_count": int(self.row_count),
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_evidence(
        cls,
        evidence: pd.DataFrame,
        *,
        splits: DonorSplitManifest,
        evidence_path: Path | str,
        source_run_manifest: RunManifest,
        feature_schema_version: str,
        purpose: str = DEVELOPMENT_EVIDENCE_PURPOSE,
    ) -> DevelopmentEvidenceProvenance:
        run_content_id, stage_content_id, dataset_content_id = _validated_calibration_source(
            evidence_path,
            splits=splits,
            source_run_manifest=source_run_manifest,
        )
        return cls(
            purpose=str(purpose),
            source_run_id=splits.source_run_id,
            split_manifest_content_id=splits.content_id,
            source_run_manifest_content_id=run_content_id,
            source_stage_manifest_content_id=stage_content_id,
            source_artifact_content_id=dataset_content_id,
            labeled_frame_fingerprint=evidence_frame_fingerprint(evidence),
            feature_schema_version=str(feature_schema_version),
            row_count=len(evidence),
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> DevelopmentEvidenceProvenance:
        if not isinstance(value, Mapping):
            raise TypeError("evidence provenance must be a mapping")
        allowed = {
            "schema_version",
            "purpose",
            "source_run_id",
            "split_manifest_content_id",
            "source_run_manifest_content_id",
            "source_stage_manifest_content_id",
            "source_artifact_content_id",
            "labeled_frame_fingerprint",
            "feature_schema_version",
            "row_count",
            "content_id",
        }
        unknown = sorted(set(value) - allowed)
        if unknown:
            raise ValueError(f"evidence provenance has unknown fields: {unknown}")
        if not str(value.get("content_id", "")).strip():
            raise ValueError("evidence provenance is missing content_id")
        row_count = value.get("row_count")
        if isinstance(row_count, bool) or not isinstance(row_count, int):
            raise TypeError("evidence provenance row_count must be a JSON integer")
        return cls(
            schema_version=int(value.get("schema_version", -1)),
            purpose=str(value.get("purpose", "")),
            source_run_id=str(value.get("source_run_id", "")),
            split_manifest_content_id=str(value.get("split_manifest_content_id", "")),
            source_run_manifest_content_id=str(value.get("source_run_manifest_content_id", "")),
            source_stage_manifest_content_id=str(value.get("source_stage_manifest_content_id", "")),
            source_artifact_content_id=str(value.get("source_artifact_content_id", "")),
            labeled_frame_fingerprint=str(value.get("labeled_frame_fingerprint", "")),
            feature_schema_version=str(value.get("feature_schema_version", "")),
            row_count=row_count,
            content_id=str(value.get("content_id", "")),
        )

    @classmethod
    def read_json(cls, path: Path | str) -> DevelopmentEvidenceProvenance:
        return cls.from_dict(_read_strict_json(path))


def _validated_calibration_source(
    evidence_path: Path | str,
    *,
    splits: DonorSplitManifest,
    source_run_manifest: RunManifest,
) -> tuple[str, str, str]:
    """Resolve an evidence dataset through its exact completed source run."""

    source_run_manifest.verify_identity()
    if source_run_manifest.run_name != splits.source_run_id:
        raise ValueError("evidence run manifest differs from the donor split source run")
    splits.validate_cohort(source_run_manifest.expected_donors)
    references = [
        reference for reference in source_run_manifest.stages if reference.stage == "calibrate"
    ]
    if len(references) != 1:
        raise ValueError("source run must contain exactly one calibrate stage reference")
    reference = references[0]
    stage = StageManifest.read_json(reference.manifest_path)
    if stage.stage != "calibrate":
        raise ValueError("source evidence manifest is not a calibrate stage")
    if (
        stage.content_id != reference.manifest_content_id
        or stage.output.content_sha256 != reference.output_content_sha256
    ):
        raise ValueError("source run calibration reference differs from its stage manifest")
    if stage.expected_donors != source_run_manifest.expected_donors:
        raise ValueError("source calibration donor set differs from its run manifest")
    source = Path(evidence_path).expanduser().resolve()
    recorded = Path(stage.output.root).expanduser().resolve()
    if source != recorded:
        raise ValueError(
            "training evidence path is not the calibration output referenced by the source run"
        )
    validate_dataset_snapshot_current(stage.output, mode="content")
    return (
        source_run_manifest.content_id,
        stage.content_id,
        stage.output.content_sha256,
    )


def validate_evidence_source(
    evidence_path: Path | str,
    *,
    provenance: DevelopmentEvidenceProvenance,
    splits: DonorSplitManifest,
    source_run_manifest: RunManifest,
) -> None:
    """Prove that fitting evidence is the source run's calibrated dataset."""

    run_id, stage_id, dataset_id = _validated_calibration_source(
        evidence_path,
        splits=splits,
        source_run_manifest=source_run_manifest,
    )
    expected = (
        provenance.source_run_manifest_content_id,
        provenance.source_stage_manifest_content_id,
        provenance.source_artifact_content_id,
    )
    if (run_id, stage_id, dataset_id) != expected:
        raise ValueError("training evidence source lineage differs from provenance")


def validate_training_evidence(
    evidence: pd.DataFrame,
    *,
    provenance: DevelopmentEvidenceProvenance,
    splits: DonorSplitManifest,
    feature_schema_version: str,
) -> None:
    if provenance.source_run_id != splits.source_run_id:
        raise ValueError("training evidence and donor split refer to different runs")
    if provenance.split_manifest_content_id != splits.content_id:
        raise ValueError("training evidence refers to a different donor split")
    if provenance.feature_schema_version != feature_schema_version:
        raise ValueError("training evidence feature schema differs from fitting")
    if provenance.row_count != len(evidence):
        raise ValueError("training evidence row count differs from provenance")
    if provenance.labeled_frame_fingerprint != evidence_frame_fingerprint(evidence):
        raise ValueError("training evidence differs from its content fingerprint")


REVIEW_LABEL_COLUMNS = (
    "donor_id",
    "object_id",
    "reference_label",
    "reviewer_1",
    "reviewer_2",
    "adjudicated",
    "root_cause",
    "evidence_sufficient",
)
REVIEW_SAMPLE_AUDIT_COLUMNS = (
    "donor_id",
    "object_id",
    "sampling_stratum",
    "sampling_method",
    "sampling_salt",
    "stratum_population",
    "stratum_sampled",
    "sampling_weight",
    "selection_probability",
)
REVIEW_SAMPLE_BLINDED_COLUMNS = (
    "donor_id",
    "object_id",
    "sampling_stratum",
)
_REVIEW_ACCEPTED_STATUSES = frozenset({"anchor", "inferred"})


def _normalise_review_target(level: object, parent: object) -> tuple[str, str]:
    target_level = str(level).strip().lower()
    target_parent = "" if parent is None else str(parent).strip()
    if target_level not in {"broad", "subtype", "specific"}:
        raise ValueError(
            "review sampling level must be 'broad', 'subtype', or 'specific'"
        )
    if (target_level in {"broad", "specific"} and target_parent) or (
        target_level == "subtype" and not target_parent
    ):
        raise ValueError(
            "broad/specific review sampling has no parent; subtype sampling requires one"
        )
    return target_level, target_parent


def _review_sampling_population(
    assignments: pd.DataFrame,
    *,
    level: str,
    parent: str | None,
    blinding_key: bytes,
    sampling_context: str,
) -> pd.DataFrame:
    """Return the exact target population and secret-keyed opaque strata."""

    if not isinstance(assignments, pd.DataFrame) or assignments.empty:
        raise ValueError("review sampling assignments must be a nonempty DataFrame")
    if assignments.columns.has_duplicates:
        raise ValueError("review sampling assignments contain duplicate columns")
    target_level, target_parent = _normalise_review_target(level, parent)
    prefix = "broad" if target_level == "broad" else "subtype"
    status_column = (
        "specific_assignment_status"
        if target_level == "specific"
        else f"{prefix}_assignment_status"
    )
    label_column = "specific_type" if target_level == "specific" else f"{prefix}_label"
    required = {
        "donor_id",
        "object_id",
        status_column,
        label_column,
    }
    if target_level == "subtype":
        required.update(
            {
                "broad_assignment_status",
                "broad_label",
                "typing_subtype_model_source",
            }
        )
    missing = sorted(required - set(assignments.columns))
    if missing:
        raise KeyError(f"review sampling assignments are missing columns: {missing}")
    work = assignments.loc[:, sorted(required)].copy()
    for column in required:
        work[column] = work[column].astype("string").fillna("").str.strip()
    if work[["donor_id", "object_id"]].eq("").any().any():
        raise ValueError("review sampling assignments contain blank donor/object keys")
    if work.duplicated(["donor_id", "object_id"]).any():
        raise ValueError("review sampling assignments contain duplicate donor/object keys")
    if target_level == "subtype":
        scope = (
            work["broad_assignment_status"].isin(_REVIEW_ACCEPTED_STATUSES)
            & work["broad_label"].eq(target_parent)
            & work["typing_subtype_model_source"].eq("probabilistic_bundle")
        )
        work = work.loc[scope].copy()
    if work.empty:
        raise ValueError(f"review sampling target {(target_level, target_parent)} is empty")
    if work[[status_column, label_column]].eq("").any().any():
        raise ValueError("review sampling target contains blank status/label strata")
    work["_raw_sampling_stratum"] = (
        work["donor_id"]
        + "|"
        + work[status_column]
        + "|"
        + work[label_column]
    )
    raw_strata = sorted(work["_raw_sampling_stratum"].unique())
    codes = {
        value: "S"
        + hmac.new(
            blinding_key,
            f"{sampling_context}|stratum|{value}".encode(),
            hashlib.sha256,
        ).hexdigest()[:20]
        for value in raw_strata
    }
    if len(set(codes.values())) != len(codes):  # pragma: no cover - cryptographic collision
        raise RuntimeError("opaque review stratum code collision")
    work["sampling_stratum"] = work["_raw_sampling_stratum"].map(codes)
    return work.loc[:, ["donor_id", "object_id", "sampling_stratum"]]


def _bound_sampling_design(
    *,
    splits: DonorSplitManifest,
    split_name: str,
    candidate_assignment_run_id: str,
    level: str,
    parent: str | None,
    blinding_key: str | bytes,
) -> tuple[str, str, str, bytes, int]:
    if not isinstance(splits, DonorSplitManifest):
        raise TypeError("splits must be DonorSplitManifest")
    split = str(split_name).strip().lower()
    if split not in REVIEW_SPLIT_NAMES:
        raise ValueError(
            "review sampling split_name must be calibration, validation, or locked_test"
        )
    assignment_run_id = str(candidate_assignment_run_id).strip().lower()
    if len(assignment_run_id) != 64 or any(
        character not in "0123456789abcdef" for character in assignment_run_id
    ):
        raise ValueError(
            "candidate_assignment_run_id must be a lowercase SHA-256 digest"
        )
    target_level, target_parent = _normalise_review_target(level, parent)
    if isinstance(blinding_key, str):
        key = blinding_key.encode("utf-8")
    elif isinstance(blinding_key, bytes):
        key = blinding_key
    else:
        raise TypeError("blinding_key must be str or bytes")
    if len(key) < 16:
        raise ValueError("blinding_key must contain at least 16 bytes")
    commitment = hashlib.sha256(key).hexdigest()
    if not splits.review_blinding_key_sha256:
        raise ValueError(
            "donor split has no precommitted review blinding-key digest"
        )
    if not hmac.compare_digest(commitment, splits.review_blinding_key_sha256):
        raise ValueError("review blinding key differs from the donor-split commitment")
    context = json.dumps(
        {
            "method": REVIEW_SAMPLING_METHOD,
            "split_manifest_content_id": splits.content_id,
            "split_name": split,
            # The assignment run identity is path/encoding independent.  Exact
            # candidate manifest/output IDs remain review-lineage fields, but
            # cannot influence which cells or bootstrap draws are selected.
            "candidate_assignment_run_id": assignment_run_id,
            "level": target_level,
            "parent": target_parent,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return target_level, target_parent, context, key, splits.review_n_per_stratum[split]


def build_review_sampling_frame(
    assignments: pd.DataFrame,
    *,
    splits: DonorSplitManifest,
    split_name: str,
    candidate_assignment_run_id: str,
    level: str,
    parent: str | None,
    blinding_key: str | bytes,
) -> pd.DataFrame:
    """Select the deterministic probability sample used for weighted review.

    The only supported probability-sample design is secret-keyed HMAC-SHA-256
    ranking within donor × assignment-status × assigned-label strata.  Sample
    size is preregistered in the donor split, and the draw is bound to the exact
    candidate manifest and target.  Use :func:`blinded_review_sample_frame`
    before handing any frame to reviewers; this returned audit frame retains
    weights and the committed sampling proof.
    """

    target_level, target_parent, context, key, requested = _bound_sampling_design(
        splits=splits,
        split_name=split_name,
        candidate_assignment_run_id=candidate_assignment_run_id,
        level=level,
        parent=parent,
        blinding_key=blinding_key,
    )
    assignment_donors = tuple(
        sorted(assignments["donor_id"].astype("string").fillna("").str.strip().unique())
    )
    expected_donors = getattr(splits, str(split_name).strip().lower())
    if assignment_donors != expected_donors:
        raise ValueError(
            "review sampling assignments must cover exactly the preregistered split donors"
        )
    population = _review_sampling_population(
        assignments,
        level=target_level,
        parent=(None if target_level == "broad" else target_parent),
        blinding_key=key,
        sampling_context=context,
    )
    population["_rank"] = [
        hmac.new(
            key,
            f"{context}|object|{donor}|{object_id}".encode(),
            hashlib.sha256,
        ).hexdigest()
        for donor, object_id in population[["donor_id", "object_id"]].itertuples(
            index=False, name=None
        )
    ]
    salt = hmac.new(
        key, f"{context}|selection".encode(), hashlib.sha256
    ).hexdigest()
    counts = population.groupby("sampling_stratum", observed=True).size().astype(int)
    selected: list[pd.DataFrame] = []
    for stratum, group in population.groupby("sampling_stratum", sort=True, observed=True):
        sampled = min(requested, len(group))
        chosen = group.sort_values(["_rank", "object_id"], kind="mergesort").head(sampled).copy()
        chosen["sampling_method"] = REVIEW_SAMPLING_METHOD
        chosen["sampling_salt"] = salt
        chosen["stratum_population"] = int(counts.loc[stratum])
        chosen["stratum_sampled"] = sampled
        chosen["sampling_weight"] = float(counts.loc[stratum]) / sampled
        chosen["selection_probability"] = sampled / float(counts.loc[stratum])
        selected.append(chosen)
    result = pd.concat(selected, ignore_index=True).drop(columns="_rank")
    return result.loc[:, list(REVIEW_SAMPLE_AUDIT_COLUMNS)].sort_values(
        ["donor_id", "object_id"], kind="mergesort", ignore_index=True
    )


def validate_review_sampling_against_assignments(
    sample_frame: pd.DataFrame,
    assignments: pd.DataFrame,
    *,
    splits: DonorSplitManifest,
    split_name: str,
    candidate_assignment_run_id: str,
    level: str,
    parent: str | None,
    blinding_key: str | bytes,
    expected_donors: Sequence[str] | None = None,
) -> None:
    """Prove review weights and selected keys against candidate output.

    By default the assignments must cover the full declared split.  Supplying
    ``expected_donors`` permits bounded donor-local replay; callers must then
    enforce complete donor coverage across their stream.
    """

    review_sample_frame_fingerprint(sample_frame)
    target_level, target_parent, context, key, requested = _bound_sampling_design(
        splits=splits,
        split_name=split_name,
        candidate_assignment_run_id=candidate_assignment_run_id,
        level=level,
        parent=parent,
        blinding_key=blinding_key,
    )
    assignment_donors = tuple(
        sorted(assignments["donor_id"].astype("string").fillna("").str.strip().unique())
    )
    split_donors = getattr(splits, str(split_name).strip().lower())
    if expected_donors is None:
        required_donors = split_donors
    else:
        required_donors = _normalise_donors(
            expected_donors, split="review validation expected_donors"
        )
        if not required_donors or not set(required_donors) <= set(split_donors):
            raise ValueError(
                "review validation expected_donors must be a nonempty subset of the split"
            )
    if assignment_donors != required_donors:
        raise ValueError(
            "review sampling assignments must cover exactly the requested split donors"
        )
    population = _review_sampling_population(
        assignments,
        level=target_level,
        parent=(None if target_level == "broad" else target_parent),
        blinding_key=key,
        sampling_context=context,
    )
    sample = sample_frame.copy()
    for column in ("donor_id", "object_id", "sampling_stratum", "sampling_salt"):
        sample[column] = sample[column].astype("string").fillna("").str.strip()
    methods = tuple(sample["sampling_method"].astype(str).unique())
    salts = tuple(sample["sampling_salt"].astype(str).unique())
    if methods != (REVIEW_SAMPLING_METHOD,) or len(salts) != 1:
        raise ValueError("review sample must use one supported sampling method and salt")
    salt = salts[0]
    expected_salt = hmac.new(
        key, f"{context}|selection".encode(), hashlib.sha256
    ).hexdigest()
    if not hmac.compare_digest(salt, expected_salt):
        raise ValueError("review sample salt differs from the precommitted sampling design")

    population_counts = population.groupby("sampling_stratum", observed=True).size().astype(int)
    sample_counts = sample.groupby("sampling_stratum", observed=True).size().astype(int)
    if set(sample_counts.index) != set(population_counts.index):
        raise ValueError("review sample does not cover every candidate assignment stratum")
    expected_sample_counts = population_counts.clip(upper=requested)
    if not sample_counts.reindex(expected_sample_counts.index).equals(expected_sample_counts):
        raise ValueError(
            "review sample counts differ from preregistered n_per_stratum"
        )
    declared_population = (
        sample.groupby("sampling_stratum", observed=True)["stratum_population"].first().astype(int)
    )
    expected_population = population_counts.reindex(declared_population.index)
    if expected_population.isna().any() or not np.array_equal(
        declared_population.to_numpy(), expected_population.to_numpy()
    ):
        raise ValueError("review sample stratum populations differ from candidate assignments")

    expected_key_to_stratum = population.set_index(["donor_id", "object_id"])[
        "sampling_stratum"
    ].to_dict()
    for donor, object_id, stratum in sample[
        ["donor_id", "object_id", "sampling_stratum"]
    ].itertuples(index=False, name=None):
        if expected_key_to_stratum.get((donor, object_id)) != stratum:
            raise ValueError("review sample key is absent from or mis-stratified for the target")

    population = population.copy()
    population["_rank"] = [
        hmac.new(
            key,
            f"{context}|object|{donor}|{object_id}".encode(),
            hashlib.sha256,
        ).hexdigest()
        for donor, object_id in population[["donor_id", "object_id"]].itertuples(
            index=False, name=None
        )
    ]
    for stratum, group in population.groupby("sampling_stratum", sort=True, observed=True):
        sampled = int(sample_counts.loc[stratum])
        expected = set(
            map(
                tuple,
                group.sort_values(["_rank", "object_id"], kind="mergesort")
                .head(sampled)[["donor_id", "object_id"]]
                .to_numpy(),
            )
        )
        observed = set(
            map(
                tuple,
                sample.loc[
                    sample["sampling_stratum"].eq(stratum), ["donor_id", "object_id"]
                ].to_numpy(),
            )
        )
        if observed != expected:
            raise ValueError("review sample is not the deterministic within-stratum draw")


def review_sample_frame_fingerprint(sample_frame: pd.DataFrame) -> str:
    """Fingerprint a prediction-blinded review sample and its sampling design."""

    if not isinstance(sample_frame, pd.DataFrame):
        raise TypeError("review sample must be a pandas DataFrame")
    if sample_frame.empty:
        raise ValueError("review sample is empty")
    if sample_frame.columns.has_duplicates:
        raise ValueError("review sample columns must be unique")
    missing = sorted(set(REVIEW_SAMPLE_AUDIT_COLUMNS) - set(sample_frame.columns))
    if missing:
        raise KeyError(f"review sample is missing sampling columns: {missing}")
    forbidden_tokens = (
        "prediction",
        "status",
        "reason",
        "candidate",
        "margin",
        "stability",
        "confidence",
        "model_score",
        "_pass",
    )
    forbidden = sorted(
        column
        for column in sample_frame.columns
        if any(token in column.lower() for token in forbidden_tokens)
        or column.lower()
        in {
            "label",
            "broad_label",
            "subtype_label",
            "specific_type",
            "cell_type",
            "broad_best_probability",
            "subtype_best_probability",
        }
        or ("probability" in column.lower() and column != "selection_probability")
    )
    if forbidden:
        raise ValueError(f"review sample is not prediction-blinded; forbidden columns: {forbidden}")
    work = sample_frame.copy()
    for column in (
        "donor_id",
        "object_id",
        "sampling_stratum",
        "sampling_method",
        "sampling_salt",
    ):
        work[column] = work[column].astype("string").fillna("").str.strip()
        if work[column].eq("").any():
            raise ValueError(f"review sample {column} contains blank values")
    if tuple(work["sampling_method"].unique()) != (REVIEW_SAMPLING_METHOD,):
        raise ValueError(
            f"review sample sampling_method must equal {REVIEW_SAMPLING_METHOD!r}"
        )
    if work["sampling_salt"].nunique() != 1:
        raise ValueError("review sample must use exactly one sampling_salt")
    if work.duplicated(["donor_id", "object_id"]).any():
        raise ValueError("review sample contains duplicate donor/object keys")
    for column in (
        "stratum_population",
        "stratum_sampled",
        "sampling_weight",
        "selection_probability",
    ):
        values = pd.to_numeric(work[column], errors="coerce")
        if (
            values.isna().any()
            or values.isin([float("inf"), float("-inf")]).any()
            or not (values > 0).all()
        ):
            raise ValueError(f"review sample {column} must be finite and positive")
        work[column] = values
    if (work["selection_probability"] > 1).any():
        raise ValueError("review sample selection_probability must lie in (0, 1]")
    population = work["stratum_population"].to_numpy(float)
    sampled = work["stratum_sampled"].to_numpy(float)
    weight = work["sampling_weight"].to_numpy(float)
    selection = work["selection_probability"].to_numpy(float)
    if (
        not np.equal(population, np.floor(population)).all()
        or not np.equal(sampled, np.floor(sampled)).all()
        or (sampled > population).any()
    ):
        raise ValueError(
            "review sample stratum_population/stratum_sampled must be integer "
            "counts with sampled <= population"
        )
    if not np.allclose(weight, population / sampled, rtol=1e-12, atol=1e-12):
        raise ValueError("review sample sampling_weight must equal population / sampled")
    if not np.allclose(selection, 1.0 / weight, rtol=1e-12, atol=1e-12):
        raise ValueError("review sample selection_probability must equal 1 / sampling_weight")
    stratum_counts = work.groupby("sampling_stratum", observed=True).size()
    design_columns = (
        "stratum_population",
        "stratum_sampled",
        "sampling_weight",
        "selection_probability",
        "sampling_method",
        "sampling_salt",
    )
    declared_counts = work.groupby("sampling_stratum", observed=True)[
        list(design_columns)
    ].nunique()
    if not declared_counts.eq(1).all().all():
        raise ValueError("review sample strata have inconsistent sampling design values")
    declared = work.groupby("sampling_stratum", observed=True)["stratum_sampled"].first()
    if not stratum_counts.astype(float).eq(declared.astype(float)).all():
        raise ValueError("review sample does not contain every selected stratum row")
    ordered = [
        "donor_id",
        "object_id",
        *sorted(set(work.columns) - {"donor_id", "object_id"}),
    ]
    return dataframe_content_sha256(work.loc[:, ordered], sort_by=["donor_id", "object_id"])


def blinded_review_sample_frame(
    sampling_audit: pd.DataFrame,
    *,
    locations: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Create the only frame that should be handed to cell reviewers.

    Sampling weights, population sizes, the committed salt, and every model
    output remain in the separately retained audit frame.  Optional location
    columns are joined by exact cell key after screening for prediction leaks.
    """

    review_sample_frame_fingerprint(sampling_audit)
    blinded = sampling_audit.loc[:, list(REVIEW_SAMPLE_BLINDED_COLUMNS)].copy()
    if locations is not None:
        if not isinstance(locations, pd.DataFrame):
            raise TypeError("locations must be a pandas DataFrame or None")
        if locations.columns.has_duplicates:
            raise ValueError("review locations contain duplicate columns")
        missing = sorted({"donor_id", "object_id"} - set(locations.columns))
        if missing:
            raise KeyError(f"review locations are missing key columns: {missing}")
        forbidden_tokens = (
            "prediction",
            "status",
            "reason",
            "candidate",
            "probability",
            "margin",
            "stability",
            "confidence",
            "model_score",
            "_pass",
            "sampling_weight",
            "stratum_population",
            "stratum_sampled",
            "sampling_salt",
            "selection_probability",
        )
        forbidden = sorted(
            column
            for column in locations.columns
            if column not in {"donor_id", "object_id"}
            and (
                any(token in column.lower() for token in forbidden_tokens)
                or column.lower()
                in {
                    "label",
                    "broad_label",
                    "subtype_label",
                    "specific_type",
                    "cell_type",
                }
            )
        )
        if forbidden:
            raise ValueError(
                "review locations contain prediction or sampling-design fields: "
                f"{forbidden}"
            )
        location_work = locations.copy()
        for column in ("donor_id", "object_id"):
            location_work[column] = (
                location_work[column].astype("string").fillna("").str.strip()
            )
        if location_work.duplicated(["donor_id", "object_id"]).any():
            raise ValueError("review locations contain duplicate cell keys")
        blinded = blinded.merge(
            location_work,
            on=["donor_id", "object_id"],
            how="left",
            validate="1:1",
            indicator=True,
        )
        if not blinded["_merge"].eq("both").all():
            raise ValueError("review locations do not cover every sampled cell")
        blinded = blinded.drop(columns="_merge")
    result = blinded.sort_values(
        ["donor_id", "object_id"], kind="mergesort", ignore_index=True
    )
    reviewer_frame_fingerprint(result)
    return result


def reviewer_frame_fingerprint(reviewer_frame: pd.DataFrame) -> str:
    """Fingerprint the exact prediction-blinded context shown to reviewers."""

    if not isinstance(reviewer_frame, pd.DataFrame) or reviewer_frame.empty:
        raise ValueError("reviewer frame must be a nonempty DataFrame")
    if reviewer_frame.columns.has_duplicates:
        raise ValueError("reviewer frame columns must be unique")
    missing = sorted(set(REVIEW_SAMPLE_BLINDED_COLUMNS) - set(reviewer_frame.columns))
    if missing:
        raise KeyError(f"reviewer frame is missing columns: {missing}")
    forbidden_tokens = (
        "prediction",
        "status",
        "reason",
        "candidate",
        "probability",
        "margin",
        "stability",
        "confidence",
        "model_score",
        "_pass",
        "sampling_weight",
        "stratum_population",
        "stratum_sampled",
        "sampling_salt",
        "sampling_method",
        "selection_probability",
    )
    forbidden = sorted(
        column
        for column in reviewer_frame.columns
        if column not in REVIEW_SAMPLE_BLINDED_COLUMNS
        and (
            any(token in column.lower() for token in forbidden_tokens)
            or column.lower()
            in {
                "label",
                "broad_label",
                "subtype_label",
                "specific_type",
                "cell_type",
            }
        )
    )
    if forbidden:
        raise ValueError(
            f"reviewer frame contains prediction or sampling-design fields: {forbidden}"
        )
    work = reviewer_frame.copy()
    for column in REVIEW_SAMPLE_BLINDED_COLUMNS:
        work[column] = work[column].astype("string").fillna("").str.strip()
        if work[column].eq("").any():
            raise ValueError(f"reviewer frame {column} contains blanks")
    if work.duplicated(["donor_id", "object_id"]).any():
        raise ValueError("reviewer frame contains duplicate donor/object keys")
    ordered = [
        "donor_id",
        "object_id",
        *sorted(set(work.columns) - {"donor_id", "object_id"}),
    ]
    return dataframe_content_sha256(
        work.loc[:, ordered], sort_by=["donor_id", "object_id"]
    )


def reviewer_frame_context_artifact_content_id(
    reviewer_frame: pd.DataFrame,
) -> str:
    """Return the immutable ingest artifact stamped on a reviewer frame.

    The context identifier is never accepted as a parallel free-text argument:
    doing so would allow a caller to fingerprint arbitrary coordinates while
    merely *claiming* that they came from the bound ingest artifact.
    """

    if not isinstance(reviewer_frame, pd.DataFrame):
        raise TypeError("reviewer_frame must be a pandas DataFrame")
    content_id = str(
        reviewer_frame.attrs.get("review_context_artifact_content_id", "")
    ).strip().lower()
    if len(content_id) != 64 or any(
        character not in "0123456789abcdef" for character in content_id
    ):
        raise ValueError(
            "reviewer frame must carry a lowercase SHA-256 "
            "review_context_artifact_content_id attribute"
        )
    return content_id


def blinded_review_sample_from_dataset(
    sampling_audit: pd.DataFrame,
    *,
    source_run_manifest: RunManifest,
    splits: DonorSplitManifest,
    location_columns: Sequence[str] = (
        "image",
        "X_centroid",
        "Y_centroid",
        "cell_region",
    ),
) -> pd.DataFrame:
    """Join reviewer locations from one content-validated source dataset.

    Only the requested columns and sampled keys are materialized, one donor at
    a time.  The snapshot is content-validated before and after the read so a
    reviewer cannot silently receive coordinates or image identifiers from an
    unrelated or concurrently changed ingest artifact.
    """

    if not isinstance(source_run_manifest, RunManifest):
        raise TypeError("source_run_manifest must be RunManifest")
    if not isinstance(splits, DonorSplitManifest):
        raise TypeError("splits must be DonorSplitManifest")
    source_run_manifest.verify_identity()
    if source_run_manifest.run_name != splits.source_run_id:
        raise ValueError("location source run differs from the donor split")
    if source_run_manifest.expected_donors != tuple(sorted(splits.donors)):
        raise ValueError("location source run donor universe differs from the donor split")
    ingest_references = [
        reference for reference in source_run_manifest.stages if reference.stage == "ingest"
    ]
    if len(ingest_references) != 1:
        raise ValueError("source run must contain exactly one ingest stage")
    ingest_manifest = ingest_references[0].validate_current(
        validate_stage=True, validation_mode="content"
    )
    if ingest_manifest.stage != "ingest":  # defensive against a malformed reference
        raise ValueError("source run ingest reference resolved to a different stage")
    location_snapshot: DatasetSnapshot = ingest_manifest.output
    if isinstance(location_columns, (str, bytes)) or not isinstance(
        location_columns, Sequence
    ):
        raise TypeError("location_columns must be a sequence")
    requested = tuple(str(column).strip() for column in location_columns)
    if any(not column for column in requested):
        raise ValueError("location_columns must contain only nonempty names")
    if len(set(requested)) != len(requested):
        raise ValueError("location_columns contains duplicates")
    if set(requested) & {"donor_id", "object_id", "sampling_stratum"}:
        raise ValueError("location_columns must not repeat review key/stratum columns")
    schema_names = {field.name for field in location_snapshot.schema}
    missing = sorted(set(("object_id", *requested)) - schema_names)
    if missing:
        raise KeyError(f"location snapshot is missing columns: {missing}")
    review_sample_frame_fingerprint(sampling_audit)
    validate_dataset_snapshot_current(location_snapshot, mode="content")
    keys = sampling_audit.loc[:, ["donor_id", "object_id"]].copy()
    for column in ("donor_id", "object_id"):
        keys[column] = keys[column].astype("string").fillna("").str.strip()
    sampled_donors = tuple(sorted(keys["donor_id"].unique()))
    snapshots = {donor.donor_id: donor for donor in location_snapshot.donors}
    unknown = sorted(set(sampled_donors) - set(snapshots))
    if unknown:
        raise ValueError(f"sample donors are absent from the location snapshot: {unknown}")
    source_has_donor = "donor_id" in schema_names
    read_columns = [
        *(("donor_id",) if source_has_donor else ()),
        "object_id",
        *requested,
    ]
    selected_frames: list[pd.DataFrame] = []
    for donor_id in sampled_donors:
        donor_snapshot = snapshots[donor_id]
        frames = [
            pd.read_parquet(fingerprint.path, columns=read_columns)
            for fingerprint in donor_snapshot.files
        ]
        if not frames:
            raise ValueError(f"location snapshot donor {donor_id!r} has no files")
        source = pd.concat(frames, ignore_index=True)
        if not source_has_donor:
            source.insert(0, "donor_id", donor_id)
        for column in ("donor_id", "object_id"):
            source[column] = source[column].astype("string").fillna("").str.strip()
        if not source["donor_id"].eq(donor_id).all():
            raise ValueError("location rows escape their donor snapshot partition")
        if source.duplicated(["donor_id", "object_id"]).any():
            raise ValueError("location snapshot contains duplicate donor/object keys")
        donor_keys = keys.loc[keys["donor_id"].eq(donor_id)]
        selected = donor_keys.merge(
            source.loc[:, ["donor_id", "object_id", *requested]],
            on=["donor_id", "object_id"],
            how="left",
            validate="1:1",
            indicator=True,
        )
        if not selected["_merge"].eq("both").all():
            raise ValueError(
                f"location snapshot does not cover every sampled cell for donor {donor_id!r}"
            )
        selected_frames.append(selected.drop(columns="_merge"))
    locations = pd.concat(selected_frames, ignore_index=True)
    result = blinded_review_sample_frame(sampling_audit, locations=locations)
    validate_dataset_snapshot_current(location_snapshot, mode="content")
    result.attrs["review_context_artifact_content_id"] = location_snapshot.content_sha256
    return result


def _canonical_review_frame(labels: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(labels, pd.DataFrame):
        raise TypeError("review labels must be a pandas DataFrame")
    missing = sorted(set(REVIEW_LABEL_COLUMNS) - set(labels.columns))
    if missing:
        raise KeyError(f"review labels are missing columns: {missing}")
    work = labels.copy()
    text_columns = (
        "donor_id",
        "object_id",
        "reference_label",
        "reviewer_1",
        "reviewer_2",
        "root_cause",
    )
    for column in text_columns:
        work[column] = work[column].astype("string").fillna("").str.strip()
    if (work[["donor_id", "object_id"]] == "").any().any():
        raise ValueError("review donor_id and object_id must be non-empty")
    if work.duplicated(["donor_id", "object_id"]).any():
        raise ValueError("review labels contain duplicate donor/object keys")
    for column in ("adjudicated", "evidence_sufficient"):
        values = work[column]
        if values.isna().any() or not values.isin([True, False, 0, 1]).all():
            raise ValueError(f"{column} must contain explicit booleans")
        work[column] = values.astype(bool)
    adjudicated = work["adjudicated"]
    required_text = (
        work["reference_label"].ne("")
        & work["reviewer_1"].ne("")
        & work["reviewer_2"].ne("")
        & work["root_cause"].ne("")
    )
    if (adjudicated & ~required_text).any():
        raise ValueError("adjudicated review rows require reference, two reviewers, and root cause")
    if (adjudicated & work["reviewer_1"].eq(work["reviewer_2"])).any():
        raise ValueError("adjudicated review rows require two distinct reviewers")
    sufficient = work["evidence_sufficient"]
    unresolved = (
        work["reference_label"].str.lower().isin({"", "unsure", "ambiguous", "unavailable"})
    )
    if (sufficient & (~adjudicated | unresolved)).any():
        raise ValueError("evidence-sufficient review rows require an adjudicated biological label")
    return work


def _review_is_complete_evidence_sufficient(labels: pd.DataFrame) -> bool:
    """Return whether every sampled row is resolved and fit for evaluation."""

    unresolved = (
        labels["reference_label"].str.lower().isin({"", "unsure", "ambiguous", "unavailable"})
    )
    return bool(
        labels["adjudicated"].all() and labels["evidence_sufficient"].all() and (~unresolved).all()
    )


def review_label_frame_fingerprint(labels: pd.DataFrame) -> str:
    """Fingerprint the complete blinded/adjudicated review ledger."""

    work = _canonical_review_frame(labels)
    records = (
        work.loc[:, list(REVIEW_LABEL_COLUMNS)]
        .sort_values(["donor_id", "object_id"], kind="mergesort")
        .to_dict(orient="records")
    )
    return sha256_json(
        {"schema_version": REVIEW_LABEL_SCHEMA_VERSION, "review_labels": records}
    )


@dataclass(frozen=True)
class ReviewLabelProvenance:
    """Manifest for one split-restricted, content-addressed review ledger."""

    bundle_id: str
    source_run_id: str
    split_manifest_content_id: str
    assignment_run_id: str
    typing_model_bundle_content_id: str
    candidate_evaluation_manifest_content_id: str
    candidate_assignment_output_content_id: str
    review_context_artifact_content_id: str
    sample_frame_fingerprint: str
    reviewer_frame_fingerprint: str
    sample_row_count: int
    split_name: str
    fingerprint: str
    purpose: str
    level: str = "broad"
    parent: str | None = None
    complete_evidence_sufficient: bool = False
    schema_version: int = REVIEW_LABEL_SCHEMA_VERSION
    content_id: str = ""

    def __post_init__(self) -> None:
        if int(self.schema_version) != REVIEW_LABEL_SCHEMA_VERSION:
            raise ValueError(f"unsupported review-label schema_version {self.schema_version!r}")
        if self.split_name not in REVIEW_LABEL_PURPOSES:
            raise ValueError(f"unknown review-label split {self.split_name!r}")
        expected_purpose = REVIEW_LABEL_PURPOSES[self.split_name]
        if self.purpose != expected_purpose:
            raise ValueError(f"{self.split_name} review purpose must be {expected_purpose!r}")
        level = str(self.level).strip().lower()
        parent = "" if self.parent is None else str(self.parent).strip()
        if level not in {"broad", "subtype", "specific"}:
            raise ValueError("review level must be 'broad', 'subtype', or 'specific'")
        if (level in {"broad", "specific"} and parent) or (
            level == "subtype" and not parent
        ):
            raise ValueError(
                "broad/specific reviews require blank parent; subtype reviews require parent"
            )
        if not isinstance(self.complete_evidence_sufficient, bool):
            raise TypeError("complete_evidence_sufficient must be a boolean")
        object.__setattr__(self, "level", level)
        object.__setattr__(self, "parent", parent)
        for field_name in (
            "bundle_id",
            "source_run_id",
            "split_manifest_content_id",
            "assignment_run_id",
            "typing_model_bundle_content_id",
            "candidate_evaluation_manifest_content_id",
            "candidate_assignment_output_content_id",
            "review_context_artifact_content_id",
            "sample_frame_fingerprint",
            "reviewer_frame_fingerprint",
            "fingerprint",
        ):
            if not str(getattr(self, field_name)).strip():
                raise ValueError(f"{field_name} must be a non-empty string")
        if (
            len(self.split_manifest_content_id) != 64
            or len(self.typing_model_bundle_content_id) != 64
            or len(self.candidate_evaluation_manifest_content_id) != 64
            or len(self.candidate_assignment_output_content_id) != 64
            or len(self.review_context_artifact_content_id) != 64
            or len(self.sample_frame_fingerprint) != 64
            or len(self.reviewer_frame_fingerprint) != 64
            or len(self.fingerprint) != 64
        ):
            raise ValueError("review split/label fingerprints must be SHA-256")
        if isinstance(self.sample_row_count, bool) or int(self.sample_row_count) < 1:
            raise ValueError("review sample_row_count must be a positive integer")
        object.__setattr__(self, "sample_row_count", int(self.sample_row_count))
        expected = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected:
            raise ValueError("review-label provenance content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": int(self.schema_version),
            "bundle_id": self.bundle_id,
            "source_run_id": self.source_run_id,
            "split_manifest_content_id": self.split_manifest_content_id,
            "assignment_run_id": self.assignment_run_id,
            "typing_model_bundle_content_id": self.typing_model_bundle_content_id,
            "candidate_evaluation_manifest_content_id": (
                self.candidate_evaluation_manifest_content_id
            ),
            "candidate_assignment_output_content_id": (
                self.candidate_assignment_output_content_id
            ),
            "review_context_artifact_content_id": (
                self.review_context_artifact_content_id
            ),
            "sample_frame_fingerprint": self.sample_frame_fingerprint,
            "reviewer_frame_fingerprint": self.reviewer_frame_fingerprint,
            "sample_row_count": self.sample_row_count,
            "split_name": self.split_name,
            "level": self.level,
            "parent": self.parent,
            "complete_evidence_sufficient": self.complete_evidence_sufficient,
            "purpose": self.purpose,
            "fingerprint": self.fingerprint,
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_labels(
        cls,
        labels: pd.DataFrame,
        *,
        bundle_id: str,
        splits: DonorSplitManifest,
        split_name: str,
        sample_frame: pd.DataFrame,
        reviewer_frame: pd.DataFrame,
        assignment_run_id: str,
        typing_model_bundle_content_id: str,
        candidate_evaluation_manifest_content_id: str,
        candidate_assignment_output_content_id: str,
        level: str = "broad",
        parent: str | None = None,
    ) -> ReviewLabelProvenance:
        split = str(split_name)
        if split not in REVIEW_LABEL_PURPOSES:
            raise ValueError(f"unknown review-label split {split!r}")
        canonical_labels = _canonical_review_frame(labels)
        sample_fingerprint = review_sample_frame_fingerprint(sample_frame)
        shown_frame_fingerprint = reviewer_frame_fingerprint(reviewer_frame)
        review_context_artifact_content_id = (
            reviewer_frame_context_artifact_content_id(reviewer_frame)
        )
        allowed = set(getattr(splits, split))
        unknown = sorted(set(canonical_labels["donor_id"]) - allowed)
        if unknown:
            raise ValueError(f"review labels escape declared {split} donors: {unknown}")
        label_keys = set(map(tuple, canonical_labels[["donor_id", "object_id"]].to_numpy()))
        sample_keys = set(
            map(
                tuple,
                sample_frame[["donor_id", "object_id"]].astype(str).to_numpy(),
            )
        )
        if label_keys != sample_keys:
            raise ValueError("review labels must cover the exact blinded sample")
        reviewer_keys = set(
            map(
                tuple,
                reviewer_frame[["donor_id", "object_id"]].astype(str).to_numpy(),
            )
        )
        if reviewer_keys != sample_keys:
            raise ValueError("reviewer frame must cover the exact sampling audit")
        return cls(
            bundle_id=str(bundle_id),
            source_run_id=splits.source_run_id,
            split_manifest_content_id=splits.content_id,
            assignment_run_id=str(assignment_run_id),
            typing_model_bundle_content_id=str(typing_model_bundle_content_id),
            candidate_evaluation_manifest_content_id=str(
                candidate_evaluation_manifest_content_id
            ),
            candidate_assignment_output_content_id=str(
                candidate_assignment_output_content_id
            ),
            review_context_artifact_content_id=str(
                review_context_artifact_content_id
            ),
            sample_frame_fingerprint=sample_fingerprint,
            reviewer_frame_fingerprint=shown_frame_fingerprint,
            sample_row_count=len(sample_frame),
            split_name=split,
            level=level,
            parent=parent,
            complete_evidence_sufficient=_review_is_complete_evidence_sufficient(canonical_labels),
            purpose=REVIEW_LABEL_PURPOSES[split],
            fingerprint=review_label_frame_fingerprint(labels),
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> ReviewLabelProvenance:
        if not isinstance(value, Mapping):
            raise TypeError("review-label provenance must be a mapping")
        allowed = {
            "schema_version",
            "bundle_id",
            "source_run_id",
            "split_manifest_content_id",
            "assignment_run_id",
            "typing_model_bundle_content_id",
            "candidate_evaluation_manifest_content_id",
            "candidate_assignment_output_content_id",
            "review_context_artifact_content_id",
            "sample_frame_fingerprint",
            "reviewer_frame_fingerprint",
            "sample_row_count",
            "split_name",
            "level",
            "parent",
            "complete_evidence_sufficient",
            "purpose",
            "fingerprint",
            "content_id",
        }
        unknown = sorted(set(value) - allowed)
        if unknown:
            raise ValueError(f"review-label provenance has unknown fields: {unknown}")
        if not str(value.get("content_id", "")).strip():
            raise ValueError("review-label provenance is missing content_id")
        required = {
            "level",
            "parent",
            "complete_evidence_sufficient",
            "candidate_evaluation_manifest_content_id",
            "candidate_assignment_output_content_id",
            "review_context_artifact_content_id",
            "reviewer_frame_fingerprint",
        }
        missing = sorted(required - set(value))
        if missing:
            raise ValueError(f"review-label provenance is missing required fields: {missing}")
        if not isinstance(value["level"], str):
            raise TypeError("review-label provenance level must be a string")
        if not isinstance(value["parent"], str):
            raise TypeError("review-label provenance parent must be a string")
        if not isinstance(value["complete_evidence_sufficient"], bool):
            raise TypeError(
                "review-label provenance complete_evidence_sufficient must be a boolean"
            )
        return cls(
            schema_version=int(value.get("schema_version", -1)),
            bundle_id=str(value.get("bundle_id", "")),
            source_run_id=str(value.get("source_run_id", "")),
            split_manifest_content_id=str(value.get("split_manifest_content_id", "")),
            assignment_run_id=str(value.get("assignment_run_id", "")),
            typing_model_bundle_content_id=str(value.get("typing_model_bundle_content_id", "")),
            candidate_evaluation_manifest_content_id=str(
                value.get("candidate_evaluation_manifest_content_id", "")
            ),
            candidate_assignment_output_content_id=str(
                value.get("candidate_assignment_output_content_id", "")
            ),
            review_context_artifact_content_id=str(
                value.get("review_context_artifact_content_id", "")
            ),
            sample_frame_fingerprint=str(value.get("sample_frame_fingerprint", "")),
            reviewer_frame_fingerprint=str(
                value.get("reviewer_frame_fingerprint", "")
            ),
            sample_row_count=int(value.get("sample_row_count", 0)),
            split_name=str(value.get("split_name", "")),
            level=value["level"],
            parent=value["parent"],
            complete_evidence_sufficient=value["complete_evidence_sufficient"],
            purpose=str(value.get("purpose", "")),
            fingerprint=str(value.get("fingerprint", "")),
            content_id=str(value.get("content_id", "")),
        )

    @classmethod
    def read_json(cls, path: Path | str) -> ReviewLabelProvenance:
        return cls.from_dict(_read_strict_json(path))


def validate_review_labels(
    labels: pd.DataFrame,
    *,
    provenance: ReviewLabelProvenance,
    splits: DonorSplitManifest,
    sample_frame: pd.DataFrame,
    reviewer_frame: pd.DataFrame,
    assignment_run_id: str,
    typing_model_bundle_content_id: str,
    candidate_evaluation_manifest_content_id: str,
    candidate_assignment_output_content_id: str,
    allow_locked_test: bool = False,
) -> pd.DataFrame:
    """Validate review identity and confine every row to its declared split."""

    work = _canonical_review_frame(labels)
    if provenance.source_run_id != splits.source_run_id:
        raise ValueError("review labels and split manifest refer to different runs")
    if provenance.split_manifest_content_id != splits.content_id:
        raise ValueError("review labels refer to a different donor split manifest")
    if provenance.assignment_run_id != str(assignment_run_id):
        raise ValueError("review labels refer to a different assignment run")
    if provenance.typing_model_bundle_content_id != str(typing_model_bundle_content_id):
        raise ValueError("review labels refer to a different typing model bundle")
    if provenance.candidate_evaluation_manifest_content_id != str(
        candidate_evaluation_manifest_content_id
    ):
        raise ValueError("review labels refer to a different candidate evaluation manifest")
    if provenance.candidate_assignment_output_content_id != str(
        candidate_assignment_output_content_id
    ):
        raise ValueError("review labels refer to a different candidate assignment output")
    if provenance.review_context_artifact_content_id != (
        reviewer_frame_context_artifact_content_id(reviewer_frame)
    ):
        raise ValueError("review labels refer to a different review-context artifact")
    if provenance.sample_row_count != len(sample_frame):
        raise ValueError("review sample row count differs from provenance")
    if provenance.sample_frame_fingerprint != review_sample_frame_fingerprint(sample_frame):
        raise ValueError("review sample differs from its blinded content fingerprint")
    if provenance.reviewer_frame_fingerprint != reviewer_frame_fingerprint(
        reviewer_frame
    ):
        raise ValueError("reviewer frame differs from its content fingerprint")
    label_keys = set(map(tuple, work[["donor_id", "object_id"]].to_numpy()))
    sample_keys = set(map(tuple, sample_frame[["donor_id", "object_id"]].astype(str).to_numpy()))
    if label_keys != sample_keys:
        raise ValueError("review labels do not cover the exact blinded sample")
    reviewer_keys = set(
        map(
            tuple,
            reviewer_frame[["donor_id", "object_id"]].astype(str).to_numpy(),
        )
    )
    if reviewer_keys != sample_keys:
        raise ValueError("reviewer frame does not cover the exact blinded sample")
    if provenance.split_name == "locked_test" and not allow_locked_test:
        raise ValueError("locked-test labels cannot be opened in this workflow")
    allowed = set(getattr(splits, provenance.split_name))
    unknown = sorted(set(work["donor_id"]) - allowed)
    if unknown:
        raise ValueError(f"review labels escape declared {provenance.split_name} donors: {unknown}")
    observed_complete = _review_is_complete_evidence_sufficient(work)
    if provenance.complete_evidence_sufficient != observed_complete:
        raise ValueError("review completeness/evidence-sufficiency flag differs from the ledger")
    if review_label_frame_fingerprint(work) != provenance.fingerprint:
        raise ValueError("review labels differ from their content fingerprint")
    return work


def validate_development_labels(
    labels: pd.DataFrame,
    *,
    provenance: DevelopmentLabelProvenance,
    splits: DonorSplitManifest,
) -> pd.DataFrame:
    """Validate label identity and forbid validation/test leakage into fitting."""

    work = _canonical_label_frame(labels)
    if label_frame_fingerprint(work) != provenance.fingerprint:
        raise ValueError("development labels differ from their fingerprint")
    if provenance.source_run_id != splits.source_run_id:
        raise ValueError("label bundle and donor split refer to different source runs")
    if provenance.split_manifest_content_id != splits.content_id:
        raise ValueError("development labels refer to a different donor split manifest")
    unknown = sorted(set(work["donor_id"]) - set(splits.donors))
    if unknown:
        raise ValueError(f"development labels contain donors outside the split: {unknown}")
    forbidden = sorted(set(work["donor_id"]) & (set(splits.validation) | set(splits.locked_test)))
    if forbidden:
        raise ValueError(f"validation/locked-test labels cannot enter model fitting: {forbidden}")
    allowed = set(splits.development) | set(splits.calibration)
    outside = sorted(set(work["donor_id"]) - allowed)
    if outside:
        raise ValueError(f"labels are outside development/calibration splits: {outside}")
    return work


__all__ = [
    "DEVELOPMENT_EVIDENCE_PURPOSE",
    "DEVELOPMENT_LABEL_PURPOSE",
    "DEVELOPMENT_LABEL_SOURCE",
    "DEFAULT_REVIEW_N_PER_STRATUM",
    "DevelopmentEvidenceProvenance",
    "DevelopmentLabelProvenance",
    "DonorSplitManifest",
    "LABEL_COLUMNS",
    "LABEL_GOVERNANCE_COLUMNS",
    "REVIEW_LABEL_COLUMNS",
    "REVIEW_LABEL_PURPOSES",
    "REVIEW_LABEL_SCHEMA_VERSION",
    "REVIEW_SAMPLE_AUDIT_COLUMNS",
    "REVIEW_SAMPLE_BLINDED_COLUMNS",
    "REVIEW_SAMPLING_METHOD",
    "REVIEW_SPLIT_NAMES",
    "ReviewLabelProvenance",
    "SPLIT_NAMES",
    "build_review_sampling_frame",
    "blinded_review_sample_frame",
    "blinded_review_sample_from_dataset",
    "evidence_frame_fingerprint",
    "label_frame_fingerprint",
    "review_label_frame_fingerprint",
    "review_sample_frame_fingerprint",
    "reviewer_frame_context_artifact_content_id",
    "reviewer_frame_fingerprint",
    "validate_review_sampling_against_assignments",
    "validate_development_labels",
    "validate_evidence_source",
    "validate_review_labels",
    "validate_training_evidence",
]
