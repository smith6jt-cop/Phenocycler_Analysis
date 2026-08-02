"""Registry-driven hierarchical cell typing from marker evidence.

This module is the single auditable implementation of identity assignment:

* broad compartments are evaluated simultaneously rather than by a residual
  priority tree;
* only identity markers may contribute to identity, while process/state markers
  are ignored and ``SST`` is reserved for an epithelial subtype;
* curated anchor, support, and exclusion directions constrain an additive
  multinomial classifier;
* seed-cell weights are balanced first by class and then by donor;
* class priors are fixed uniform (there is no fitted intercept);
* stability is estimated by resampling whole donors, never individual cells;
* subtypes are evaluated only after a parent broad type has been accepted.

The canonical input is one row per ``donor_id``/``object_id``/``marker``.
Production evidence uses a categorical state plus ``log2_threshold_ratio``.
Typing converts that signed threshold coordinate to bounded nonnegative
support and masks invalid/indeterminate evidence before computing logits.
Empirical background-tail p-values and their ``1 - p`` summaries remain audit
quantities, not expression posteriors.  Explicit ``evidence_probability``
without a calibrated state/threshold coordinate remains supported for small
numeric-only pilots.  The output retains a best candidate for every cell with
usable identity evidence even when the authoritative assignment is
``Ambiguous``, ``Other``, or ``Unavailable``.

Registry schema
---------------
The public APIs accept either :class:`TypingRegistry` or a mapping like::

    {
        "marker_roles": {
            "CD3e": "identity",
            "Ki67": "process",
            "SST": "subtype_only",
        },
        "broad": {
            "Immune": {
                "anchors": {"CD3e": 4.0},
                "support": {"CD68": 2.0},
                "exclusions": {"E_cadherin": 2.0},
                "required": ["CD3e", "CD68"],
            },
            "Epithelial": {
                "anchors": ["E_cadherin"],
                "support": ["Pan_Cytokeratin"],
                "exclusions": ["CD3e"],
                "required": ["E_cadherin"],
            },
        },
        "subtypes": {
            "T_cell": {"parent": "Immune", "anchors": ["CD3e"]},
            "Delta": {"parent": "Epithelial", "anchors": ["SST"]},
        },
    }

Numeric rule values are transparent starting/default weights.  A sequence uses
4 for anchors and 2 for support/exclusion.  Fitting may change their magnitudes,
but never their curated signs, and unrelated marker/class coefficients remain
exactly zero.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logsumexp

DONOR_COLUMN = "donor_id"
OBJECT_COLUMN = "object_id"
MARKER_COLUMN = "marker"
KEY_COLUMNS = (DONOR_COLUMN, OBJECT_COLUMN)
DEFAULT_TYPING_RULES_PATH = Path(__file__).with_name("typing_rules.json")
FEATURE_SCHEMA_VERSION = "state-safe-threshold-support-v1"

# Export-facing spelling.  ``anchor`` means an authoritative curated anchor,
# not merely the classifier's largest coefficient.
AUTHORITATIVE_ANCHOR = "anchor"
INFERRED = "inferred"
AMBIGUOUS = "ambiguous"
OTHER = "other"
UNAVAILABLE = "unavailable"

OTHER_LABEL = "Other"
AMBIGUOUS_LABEL = "Ambiguous"
UNAVAILABLE_LABEL = "Unavailable"

IDENTITY_ROLE = "identity"
SUBTYPE_ONLY_ROLE = "subtype_only"
PROCESS_ROLE = "process"
EXCLUDED_ROLE = "excluded"

# This guard is deliberately stronger than a user-supplied registry: SST may
# subtype epithelial/endocrine cells, but cannot itself pull a cell into the
# broad Epithelial compartment.
SUBTYPE_ONLY_MARKERS = frozenset({"SST"})

_POSITIVE_STATES = frozenset({"positive", "pos", "expressed", "present", "1", "true"})
_NEGATIVE_STATES = frozenset({"negative", "neg", "absent", "0", "false"})
_INDETERMINATE_STATES = frozenset({"indeterminate", "ambiguous", "uncertain", "borderline"})
_UNAVAILABLE_STATES = frozenset({"unavailable", "missing", "invalid", "not_available", "na"})


def _as_string_tuple(value: Any) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,)
    return tuple(str(v) for v in value)


def _weighted_markers(value: Any, *, default: float) -> tuple[tuple[str, float], ...]:
    """Normalize a marker sequence/mapping into ordered ``(marker, weight)`` pairs."""
    if value is None:
        return ()
    items = value.items() if isinstance(value, Mapping) else ((m, default) for m in _as_string_tuple(value))
    out: list[tuple[str, float]] = []
    seen: set[str] = set()
    for marker, weight in items:
        marker = str(marker)
        weight = float(weight)
        if not marker:
            raise ValueError("marker names must be non-empty")
        if not np.isfinite(weight) or weight <= 0:
            raise ValueError(f"rule weight for {marker!r} must be finite and > 0")
        if marker in seen:
            raise ValueError(f"marker {marker!r} appears twice in one rule")
        seen.add(marker)
        out.append((marker, weight))
    return tuple(out)


def _normalise_role(role: Any) -> str:
    role = str(role or IDENTITY_ROLE).strip().lower().replace("-", "_")
    aliases = {
        "type": IDENTITY_ROLE,
        "lineage": IDENTITY_ROLE,
        "state": PROCESS_ROLE,
        "functional": PROCESS_ROLE,
        "subtype": SUBTYPE_ONLY_ROLE,
        "subtypeonly": SUBTYPE_ONLY_ROLE,
        "exclude": EXCLUDED_ROLE,
        "mask": EXCLUDED_ROLE,
    }
    role = aliases.get(role, role)
    if role not in {IDENTITY_ROLE, SUBTYPE_ONLY_ROLE, PROCESS_ROLE, EXCLUDED_ROLE}:
        raise ValueError(
            f"unknown marker role {role!r}; expected identity, subtype_only, process, or excluded"
        )
    return role


def _canonical_fingerprint(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True)
class TypingRulesOverlay:
    """Versioned biological rules layered over measurement calibration."""

    schema_version: int
    rules_version: str
    marker_registry_version: str
    description: str
    broad: Mapping[str, Mapping[str, Any]]
    subtypes: Mapping[str, Mapping[str, Any]]
    fingerprint: str


def typing_rules_from_dict(payload: Mapping[str, Any]) -> TypingRulesOverlay:
    """Parse and structurally validate a typing-rules mapping."""
    if not isinstance(payload, Mapping):
        raise TypeError("typing rules must be a mapping")
    required = {
        "schema_version",
        "rules_version",
        "marker_registry_version",
        "broad",
        "subtypes",
    }
    missing = sorted(required - set(payload))
    if missing:
        raise ValueError(f"typing rules are missing required keys: {missing}")
    if int(payload["schema_version"]) != 1:
        raise ValueError(f"unsupported typing-rules schema_version {payload['schema_version']!r}")
    rules_version = str(payload["rules_version"]).strip()
    registry_version = str(payload["marker_registry_version"]).strip()
    if not rules_version or not registry_version:
        raise ValueError("typing rules and marker registry versions must be non-empty")
    broad = payload["broad"]
    subtypes = payload["subtypes"]
    if not isinstance(broad, Mapping) or not broad:
        raise ValueError("typing rules 'broad' must be a non-empty mapping")
    if not isinstance(subtypes, Mapping):
        raise ValueError("typing rules 'subtypes' must be a mapping")
    for level, rules in (("broad", broad), ("subtype", subtypes)):
        for name, rule in rules.items():
            if not str(name).strip() or not isinstance(rule, Mapping):
                raise ValueError(f"{level} rule names must be non-empty and map to objects")
            _weighted_markers(rule.get("anchors", rule.get("anchor")), default=4.0)
            _weighted_markers(rule.get("support", rule.get("supports")), default=2.0)
            _weighted_markers(
                rule.get("exclusions", rule.get("exclusion", rule.get("exclude"))),
                default=2.0,
            )
            if level == "subtype" and not str(rule.get("parent", "")).strip():
                raise ValueError(f"subtype {name!r} must declare a parent")
    return TypingRulesOverlay(
        schema_version=1,
        rules_version=rules_version,
        marker_registry_version=registry_version,
        description=str(payload.get("description", "")),
        broad={str(name): dict(rule) for name, rule in broad.items()},
        subtypes={str(name): dict(rule) for name, rule in subtypes.items()},
        fingerprint=_canonical_fingerprint(payload),
    )


def load_typing_rules(path: str | Path | None = None) -> TypingRulesOverlay:
    """Load the default or explicitly selected biological typing overlay."""
    source = DEFAULT_TYPING_RULES_PATH if path is None else Path(path)
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"typing rules not found: {source}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON in typing rules {source}: {exc}") from exc
    return typing_rules_from_dict(payload)


@dataclass(frozen=True)
class ClassRule:
    """Curated additive evidence directions for one broad class or subtype."""

    name: str
    parent: str | None = None
    anchors: tuple[tuple[str, float], ...] = ()
    support: tuple[tuple[str, float], ...] = ()
    exclusions: tuple[tuple[str, float], ...] = ()
    required: tuple[str, ...] = ()

    @property
    def anchor_markers(self) -> tuple[str, ...]:
        return tuple(marker for marker, _ in self.anchors)

    @property
    def support_markers(self) -> tuple[str, ...]:
        return tuple(marker for marker, _ in self.support)

    @property
    def exclusion_markers(self) -> tuple[str, ...]:
        return tuple(marker for marker, _ in self.exclusions)

    @property
    def model_markers(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys(self.anchor_markers + self.support_markers + self.exclusion_markers))

    def initial_coefficient(self, marker: str) -> float:
        positive = [weight for m, weight in self.anchors + self.support if m == marker]
        negative = [weight for m, weight in self.exclusions if m == marker]
        if positive and negative:
            raise ValueError(f"{self.name}: marker {marker!r} is both supporting and excluding")
        if positive:
            return max(positive)
        if negative:
            return -max(negative)
        return 0.0


@dataclass(frozen=True)
class TypingRegistry:
    """Validated broad/subtype rules and marker roles."""

    broad_rules: tuple[ClassRule, ...]
    subtype_rules: tuple[ClassRule, ...] = ()
    marker_roles: Mapping[str, str] = field(default_factory=dict)
    rules_version: str = "inline-rules"
    rules_fingerprint: str = ""
    marker_registry_fingerprint: str = ""

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any] | "TypingRegistry") -> "TypingRegistry":
        if isinstance(value, cls):
            return value
        if not isinstance(value, Mapping):
            raise TypeError("registry must be a TypingRegistry or mapping")

        raw_roles: dict[str, Any] = dict(value.get("marker_roles", {}))
        for marker, spec in dict(value.get("markers", {})).items():
            if isinstance(spec, Mapping):
                raw_roles.setdefault(str(marker), spec.get("role", IDENTITY_ROLE))
            else:
                raw_roles.setdefault(str(marker), spec)
        roles = {str(marker): _normalise_role(role) for marker, role in raw_roles.items()}
        for marker in SUBTYPE_ONLY_MARKERS:
            roles[marker] = SUBTYPE_ONLY_ROLE

        raw_broad = value.get("broad", value.get("broad_classes"))
        if not isinstance(raw_broad, Mapping) or not raw_broad:
            raise ValueError("registry must define at least one class under 'broad'")

        broad: list[ClassRule] = []
        for name, spec in raw_broad.items():
            broad.append(_parse_rule(str(name), spec, parent=None, roles=roles, level="broad"))

        raw_subtypes = value.get("subtypes", value.get("subtype_classes", {}))
        if raw_subtypes is None:
            raw_subtypes = {}
        if not isinstance(raw_subtypes, Mapping):
            raise ValueError("registry 'subtypes' must be a mapping")
        subtype: list[ClassRule] = []
        for name, spec in raw_subtypes.items():
            if _looks_like_rule(spec):
                if not isinstance(spec, Mapping) or not spec.get("parent"):
                    raise ValueError(f"subtype {name!r} must declare a parent")
                subtype.append(
                    _parse_rule(
                        str(name),
                        spec,
                        parent=str(spec["parent"]),
                        roles=roles,
                        level="subtype",
                    )
                )
            else:
                # Also accept {"Epithelial": {"Beta": {...}, "Delta": {...}}}.
                if not isinstance(spec, Mapping):
                    raise ValueError(f"nested subtype block {name!r} must be a mapping")
                for subtype_name, subtype_spec in spec.items():
                    subtype.append(
                        _parse_rule(
                            str(subtype_name),
                            subtype_spec,
                            parent=str(name),
                            roles=roles,
                            level="subtype",
                        )
                    )

        broad_names = [rule.name for rule in broad]
        if len(set(broad_names)) != len(broad_names):
            raise ValueError("broad class names must be unique")
        subtype_names = [rule.name for rule in subtype]
        if len(set(subtype_names)) != len(subtype_names):
            raise ValueError("subtype names must be unique")
        unknown_parents = sorted({rule.parent for rule in subtype} - set(broad_names))
        if unknown_parents:
            raise ValueError(f"subtypes refer to unknown broad parents: {unknown_parents}")

        fingerprint = str(value.get("rules_fingerprint", "")).strip()
        if not fingerprint:
            fingerprint = _canonical_fingerprint(value)
        return cls(
            broad_rules=tuple(broad),
            subtype_rules=tuple(subtype),
            marker_roles=roles,
            rules_version=str(value.get("rules_version", "inline-rules")),
            rules_fingerprint=fingerprint,
            marker_registry_fingerprint=str(value.get("marker_registry_fingerprint", "")),
        )

    @classmethod
    def from_marker_registry(
        cls,
        marker_registry: Any,
        **kwargs: Any,
    ) -> "TypingRegistry":
        """Adapt :class:`marker_registry.MarkerRegistry` without copying marker lists.

        See :func:`typing_registry_from_marker_registry` for optional curated
        anchor, exclusion, and subtype overlays.
        """
        return typing_registry_from_marker_registry(marker_registry, **kwargs)

    def marker_role(self, marker: str) -> str:
        if marker in SUBTYPE_ONLY_MARKERS:
            return SUBTYPE_ONLY_ROLE
        return _normalise_role(self.marker_roles.get(marker, IDENTITY_ROLE))

    def rules_for(self, *, level: str = "broad", parent: str | None = None) -> tuple[ClassRule, ...]:
        if level == "broad":
            if parent is not None:
                raise ValueError("parent is only valid for subtype rules")
            return self.broad_rules
        if level != "subtype":
            raise ValueError("level must be 'broad' or 'subtype'")
        if parent is None:
            raise ValueError("a subtype parent is required")
        return tuple(rule for rule in self.subtype_rules if rule.parent == parent)


def _looks_like_rule(spec: Any) -> bool:
    if not isinstance(spec, Mapping):
        return False
    return bool(
        set(spec)
        & {
            "parent",
            "anchor",
            "anchors",
            "support",
            "supports",
            "exclude",
            "exclusion",
            "exclusions",
            "required",
            "required_markers",
            "required_models",
        }
    )


def _parse_rule(
    name: str,
    spec: Any,
    *,
    parent: str | None,
    roles: Mapping[str, str],
    level: str,
) -> ClassRule:
    if not isinstance(spec, Mapping):
        raise ValueError(f"class rule {name!r} must be a mapping")
    anchors = _weighted_markers(spec.get("anchors", spec.get("anchor")), default=4.0)
    support = _weighted_markers(spec.get("support", spec.get("supports")), default=2.0)
    exclusions = _weighted_markers(
        spec.get("exclusions", spec.get("exclusion", spec.get("exclude"))),
        default=2.0,
    )

    def role(marker: str) -> str:
        if marker in SUBTYPE_ONLY_MARKERS:
            return SUBTYPE_ONLY_ROLE
        return _normalise_role(roles.get(marker, IDENTITY_ROLE))

    allowed_roles = {IDENTITY_ROLE} if level == "broad" else {IDENTITY_ROLE, SUBTYPE_ONLY_ROLE}

    def allowed(weighted: tuple[tuple[str, float], ...]) -> tuple[tuple[str, float], ...]:
        return tuple((marker, weight) for marker, weight in weighted if role(marker) in allowed_roles)

    anchors = allowed(anchors)
    support = allowed(support)
    exclusions = allowed(exclusions)
    positive = {marker for marker, _ in anchors + support}
    negative = {marker for marker, _ in exclusions}
    overlap = sorted(positive & negative)
    if overlap:
        raise ValueError(f"{name}: marker(s) both support and exclude the class: {overlap}")
    if not positive and not negative:
        raise ValueError(f"{name}: no eligible identity markers remain after role filtering")

    raw_required = spec.get(
        "required",
        spec.get("required_markers", spec.get("required_models")),
    )
    if raw_required is None:
        # Default to the positive model inputs. Explicit ``required`` is preferred
        # whenever support markers are optional.
        required = tuple(dict.fromkeys(marker for marker, _ in anchors + support))
    else:
        required = tuple(
            marker for marker in _as_string_tuple(raw_required) if role(marker) in allowed_roles
        )
    if not required:
        raise ValueError(f"{name}: at least one eligible required marker is needed")
    unknown_required = sorted(set(required) - (positive | negative))
    if unknown_required:
        raise ValueError(
            f"{name}: required markers must also appear in anchors/support/exclusions: "
            f"{unknown_required}"
        )
    return ClassRule(
        name=name,
        parent=parent,
        anchors=anchors,
        support=support,
        exclusions=exclusions,
        required=tuple(dict.fromkeys(required)),
    )


def _resolve_key_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Accept concise ``donor``/``object`` aliases but standardize the public output."""
    out = frame.copy()
    aliases = {"donor": DONOR_COLUMN, "object": OBJECT_COLUMN}
    for alias, canonical in aliases.items():
        if canonical not in out and alias in out:
            out = out.rename(columns={alias: canonical})
    missing = [c for c in (*KEY_COLUMNS, MARKER_COLUMN) if c not in out]
    if missing:
        raise KeyError(f"cell-marker evidence is missing required columns {missing}")
    return out


def _canonical_states(series: pd.Series) -> pd.Series:
    # Calibration deliberately stores this low-cardinality field as a
    # categorical.  Convert to pandas' nullable string dtype before filling
    # missing values; filling a Categorical with ``""`` is illegal unless that
    # sentinel was predeclared as a category.
    raw = (
        series.astype("string")
        .fillna("")
        .str.strip()
        .str.lower()
        .str.replace(r"[\s-]+", "_", regex=True)
    )
    state = pd.Series("", index=series.index, dtype=object)
    state[raw.isin(_POSITIVE_STATES)] = "positive"
    state[raw.isin(_NEGATIVE_STATES)] = "negative"
    state[raw.isin(_INDETERMINATE_STATES)] = "indeterminate"
    state[raw.isin(_UNAVAILABLE_STATES)] = "unavailable"
    unknown = sorted(set(raw[(raw != "") & (state == "")]))
    if unknown:
        raise ValueError(f"unknown cell-marker state value(s): {unknown}")
    return state


def _threshold_relative_support(values: pd.Series) -> pd.Series:
    """Map a signed log2 threshold coordinate to bounded positive support.

    The centered logistic transform

    ``max(0, 2 * sigmoid(log(2) * ratio) - 1)``

    is zero at and below the donor-local threshold, remains monotone above the
    threshold, and approaches one without allowing extreme intensities to
    dominate an additive logit.  It intentionally is not called an expression
    posterior.
    """

    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    if np.isinf(numeric.to_numpy()).any():
        raise ValueError("log2 threshold-relative evidence must be finite")
    support = pd.Series(np.nan, index=numeric.index, dtype=float)
    finite = numeric.notna()
    transformed = 2.0 * expit(np.log(2.0) * numeric.loc[finite]) - 1.0
    support.loc[finite] = np.maximum(transformed, 0.0)
    return support


def _state_safe_marker_support(
    frame: pd.DataFrame,
    state: pd.Series,
    *,
    ratio_columns: Sequence[str],
    score_columns: Sequence[str],
    tail_columns: Sequence[str],
    audit_only_columns: Sequence[str],
    explicit_valid_column: str | None,
    tail_probability_is_pvalue: bool,
    context: str,
) -> pd.DataFrame:
    """Resolve one marker to a state-safe feature and audit diagnostics."""

    all_columns = tuple(
        dict.fromkeys(
            (*ratio_columns, *score_columns, *tail_columns, *audit_only_columns)
        )
    )
    numeric: dict[str, pd.Series] = {}
    for column in all_columns:
        if column not in frame:
            continue
        values = pd.to_numeric(frame[column], errors="coerce").astype(float)
        if np.isinf(values.to_numpy()).any():
            raise ValueError(f"{context}: {column} must be finite where present")
        numeric[column] = values

    probability_like = [
        column
        for column in (*score_columns, *tail_columns, *audit_only_columns)
        if column in numeric
    ]
    for column in probability_like:
        values = numeric[column]
        outside = values.notna() & ((values < 0.0) | (values > 1.0))
        if outside.any():
            bad = values.loc[outside].head(5).tolist()
            raise ValueError(
                f"{context}: {column} values must lie in [0, 1]; examples: {bad}"
            )

    ratio = pd.Series(np.nan, index=frame.index, dtype=float)
    source = pd.Series("", index=frame.index, dtype=object)
    for column in ratio_columns:
        if column not in numeric:
            continue
        candidate = numeric[column]
        take = ratio.isna() & candidate.notna()
        ratio.loc[take] = candidate.loc[take]
        source.loc[take] = column

    raw_support = _threshold_relative_support(ratio)
    no_ratio = ratio.isna()
    positive_fallback = no_ratio & state.eq("positive")
    negative_fallback = no_ratio & state.eq("negative")
    raw_support.loc[positive_fallback] = 1.0
    raw_support.loc[negative_fallback] = 0.0
    source.loc[(source == "") & (positive_fallback | negative_fallback)] = "state"

    # ``evidence_probability`` remains an explicit pilot interface only when a
    # calibrated threshold coordinate and categorical decision are both absent.
    # A tail p-value can enter this path solely when the caller explicitly says
    # it is already a score (``tail_probability_is_pvalue=False``).
    usable_scores = list(score_columns)
    if not tail_probability_is_pvalue:
        usable_scores.extend(tail_columns)
    for column in usable_scores:
        if column not in numeric:
            continue
        candidate = numeric[column]
        take = raw_support.isna() & state.eq("") & candidate.notna()
        raw_support.loc[take] = candidate.loc[take]
        source.loc[take] = column

    # Keep an explicit audit trail when a background-tail p-value is the only
    # numeric quantity.  It is deliberately not converted to ``1 - p``.
    if tail_probability_is_pvalue:
        for column in tail_columns:
            if column not in numeric:
                continue
            take = (source == "") & numeric[column].notna()
            source.loc[take] = f"{column}_audit_only"
    for column in audit_only_columns:
        if column not in numeric:
            continue
        take = (source == "") & numeric[column].notna()
        source.loc[take] = f"{column}_audit_only"

    # If a state blocks a numeric pilot score, retain the source for diagnosis
    # even though the value cannot enter a logit.
    for column in score_columns:
        if column not in numeric:
            continue
        take = (source == "") & numeric[column].notna()
        source.loc[take] = column
    source.loc[(source == "") & state.ne("")] = "state"

    numeric_available = pd.Series(False, index=frame.index, dtype=bool)
    for values in numeric.values():
        numeric_available |= values.notna()
    available = state.ne("unavailable") & (state.ne("") | numeric_available)

    model_valid = raw_support.notna() & ~state.isin(
        {"indeterminate", "unavailable"}
    )
    if explicit_valid_column is not None and explicit_valid_column in frame:
        explicit = frame[explicit_valid_column]
        if explicit.isna().any() or not explicit.isin([True, False, 0, 1]).all():
            raise ValueError(
                f"{context}: {explicit_valid_column} must contain explicit booleans"
            )
        model_valid &= explicit.astype(bool)

    # A categorical negative is authoritative for feature construction even if
    # a malformed pilot row also carries a large numeric score. Invalid and
    # indeterminate-but-measurable evidence is represented as zero; unavailable
    # or wholly absent evidence remains NaN for the assignment-state contract.
    support = raw_support.copy()
    support.loc[state.eq("negative")] = 0.0
    support.loc[available & ~model_valid] = 0.0
    support.loc[~available] = np.nan
    logit_eligible = model_valid & support.notna()

    return pd.DataFrame(
        {
            "typing_support": support,
            "evidence_probability": support,
            "evidence_available": available,
            "model_valid": model_valid,
            "logit_eligible": logit_eligible,
            "evidence_source": source,
        },
        index=frame.index,
    )


def prepare_cell_marker_evidence(
    evidence: pd.DataFrame,
    *,
    tail_probability_is_pvalue: bool = True,
) -> pd.DataFrame:
    """Normalize long-form evidence to state-safe typing support.

    Production precedence is ``log2_threshold_ratio`` (or its pilot alias
    ``log2_ratio``), then categorical positive/negative fallback.  An explicit
    numeric ``evidence_probability`` is used only when neither calibrated
    source is present, preserving numeric-only pilots. ``expression_probability``
    (the calibration artifact's ``1 - p`` field) and default upper-tail
    p-values are validated and retained as audit sources but never used as
    expression posteriors.

    Invalid and indeterminate evidence is zero in ``typing_support`` and
    ``evidence_probability`` so it cannot affect a logit.  Unavailable or absent
    evidence remains NaN. ``logit_eligible`` makes this pass/fail distinction
    explicit for diagnostics.
    """
    if not isinstance(evidence, pd.DataFrame):
        raise TypeError("evidence must be a pandas DataFrame")
    d = _resolve_key_columns(evidence)
    d[DONOR_COLUMN] = d[DONOR_COLUMN].astype(str)
    d[OBJECT_COLUMN] = d[OBJECT_COLUMN].astype(str)
    d[MARKER_COLUMN] = d[MARKER_COLUMN].astype(str)
    if (d[MARKER_COLUMN] == "").any():
        raise ValueError("marker names must be non-empty")
    duplicate = d.duplicated([*KEY_COLUMNS, MARKER_COLUMN], keep=False)
    if duplicate.any():
        sample = d.loc[duplicate, [*KEY_COLUMNS, MARKER_COLUMN]].head(5).to_dict("records")
        raise ValueError(f"duplicate donor/object/marker evidence rows; examples: {sample}")

    state = _canonical_states(d["state"]) if "state" in d else pd.Series("", index=d.index)
    contract = _state_safe_marker_support(
        d,
        state,
        ratio_columns=("log2_threshold_ratio", "log2_ratio"),
        score_columns=("evidence_probability",),
        tail_columns=("empirical_tail_probability", "tail_probability"),
        audit_only_columns=("expression_probability",),
        explicit_valid_column="model_valid",
        tail_probability_is_pvalue=tail_probability_is_pvalue,
        context="cell-marker evidence",
    )

    return pd.DataFrame(
        {
            DONOR_COLUMN: d[DONOR_COLUMN],
            OBJECT_COLUMN: d[OBJECT_COLUMN],
            MARKER_COLUMN: d[MARKER_COLUMN],
            "state": state,
            **{column: contract[column] for column in contract.columns},
        }
    )


def _as_registry(registry: Mapping[str, Any] | TypingRegistry) -> TypingRegistry:
    return TypingRegistry.from_mapping(registry)


def _overlay_rule_markers(
    rule: Mapping[str, Any],
) -> tuple[set[str], set[str], set[str], set[str]]:
    anchors = {
        marker
        for marker, _ in _weighted_markers(
            rule.get("anchors", rule.get("anchor")),
            default=4.0,
        )
    }
    support = {
        marker
        for marker, _ in _weighted_markers(
            rule.get("support", rule.get("supports")),
            default=2.0,
        )
    }
    exclusions = {
        marker
        for marker, _ in _weighted_markers(
            rule.get("exclusions", rule.get("exclusion", rule.get("exclude"))),
            default=2.0,
        )
    }
    required = set(
        _as_string_tuple(
            rule.get(
                "required",
                rule.get("required_markers", rule.get("required_models")),
            )
        )
    )
    return anchors, support, exclusions, required


def validate_typing_rules(
    overlay: TypingRulesOverlay,
    marker_registry: Any,
) -> None:
    """Validate biological rules against a concrete calibration registry.

    This is intentionally strict: all active broad marker models must appear
    exactly once as positive evidence in their registered compartment, process
    markers may appear nowhere, SST may appear only as Delta evidence, and the
    overlay is locked to one marker-registry version.
    """
    marker_specs = getattr(marker_registry, "markers", None)
    registry_version = str(getattr(marker_registry, "registry_version", ""))
    by_name = getattr(marker_registry, "by_name", None)
    if marker_specs is None or not isinstance(by_name, Mapping):
        raise TypeError("marker_registry must expose markers, by_name, and registry_version")
    errors: list[str] = []
    if overlay.marker_registry_version != registry_version:
        errors.append(
            f"overlay expects marker registry {overlay.marker_registry_version!r}, "
            f"got {registry_version!r}"
        )

    expected_broad: dict[str, set[str]] = {}
    for marker in marker_specs:
        if (
            str(getattr(marker, "kind", "")) == "type"
            and bool(getattr(marker, "active", False))
            and "broad" in tuple(getattr(marker, "uses", ()))
        ):
            compartment = getattr(marker, "broad_compartment", None)
            if compartment is not None:
                expected_broad.setdefault(str(compartment), set()).add(str(marker.name))

    missing_classes = sorted(set(expected_broad) - set(overlay.broad))
    extra_classes = sorted(set(overlay.broad) - set(expected_broad))
    if missing_classes:
        errors.append(f"broad rules missing calibrated compartments {missing_classes}")
    if extra_classes:
        errors.append(f"broad rules contain unknown compartments {extra_classes}")

    sst_locations: list[tuple[str, str, str]] = []

    def validate_marker(
        marker_name: str,
        *,
        location: str,
        require_subtype_use: bool,
    ) -> None:
        spec = by_name.get(marker_name)
        if spec is None:
            errors.append(f"{location}: unknown marker {marker_name!r}")
            return
        if str(getattr(spec, "kind", "")) != "type":
            errors.append(f"{location}: process/state marker {marker_name!r} cannot drive identity")
        if not bool(getattr(spec, "active", False)):
            errors.append(f"{location}: marker {marker_name!r} has no active calibration")
        if require_subtype_use and "subtype" not in tuple(getattr(spec, "uses", ())):
            errors.append(f"{location}: marker {marker_name!r} is not registered for subtype use")

    for class_name, rule in overlay.broad.items():
        anchors, support, exclusions, required = _overlay_rule_markers(rule)
        positive = anchors | support
        overlap = sorted(positive & exclusions)
        if overlap:
            errors.append(f"broad {class_name}: markers both support and exclude: {overlap}")
        expected = expected_broad.get(class_name, set())
        missing = sorted(expected - positive)
        extra = sorted(positive - expected)
        if missing:
            errors.append(f"broad {class_name}: active markers missing from anchors/support: {missing}")
        if extra:
            errors.append(
                f"broad {class_name}: positive markers belong to another use/compartment: {extra}"
            )
        if required and not required.issubset(positive | exclusions):
            errors.append(
                f"broad {class_name}: required markers are absent from the rule: "
                f"{sorted(required - positive - exclusions)}"
            )
        for marker in positive | exclusions:
            validate_marker(
                marker,
                location=f"broad {class_name}",
                require_subtype_use=False,
            )
            if marker == "SST":
                sst_locations.append(("broad", class_name, "positive" if marker in positive else "exclusion"))

    for subtype_name, rule in overlay.subtypes.items():
        parent = str(rule.get("parent", ""))
        if parent not in expected_broad:
            errors.append(f"subtype {subtype_name}: unknown parent {parent!r}")
        anchors, support, exclusions, required = _overlay_rule_markers(rule)
        positive = anchors | support
        if not positive:
            errors.append(f"subtype {subtype_name}: at least one anchor/support marker is required")
        overlap = sorted(positive & exclusions)
        if overlap:
            errors.append(f"subtype {subtype_name}: markers both support and exclude: {overlap}")
        if required and not required.issubset(positive | exclusions):
            errors.append(
                f"subtype {subtype_name}: required markers are absent from the rule: "
                f"{sorted(required - positive - exclusions)}"
            )
        for marker in positive:
            validate_marker(
                marker,
                location=f"subtype {subtype_name}",
                require_subtype_use=True,
            )
            if marker == "SST":
                sst_locations.append(("subtype", subtype_name, "positive"))
        for marker in exclusions:
            validate_marker(
                marker,
                location=f"subtype {subtype_name} exclusion",
                require_subtype_use=False,
            )
            if marker == "SST":
                sst_locations.append(("subtype", subtype_name, "exclusion"))

    if sst_locations != [("subtype", "Delta", "positive")]:
        errors.append(
            "SST must appear exactly once, as positive evidence for the Delta subtype; "
            f"found {sst_locations}"
        )
    immune = overlay.broad.get("Immune", {})
    immune_positive = set()
    if immune:
        immune_anchors, immune_support, _immune_exclusions, _immune_required = (
            _overlay_rule_markers(immune)
        )
        immune_positive = immune_anchors | immune_support
    cd11b = by_name.get("CD11b")
    if "CD11b" not in immune_positive:
        errors.append("CD11b must remain positive broad Immune evidence")
    if cd11b is None or getattr(cd11b, "reference", None) != "E_cadherin":
        errors.append("CD11b must retain E_cadherin as its calibration reference")

    if errors:
        raise ValueError("invalid typing rules:\n- " + "\n- ".join(errors))


def typing_registry_from_marker_registry(
    marker_registry: Any,
    *,
    typing_rules: TypingRulesOverlay | Mapping[str, Any] | str | Path | None = None,
    anchors: Mapping[str, Sequence[str] | Mapping[str, float]] | None = None,
    exclusions: Mapping[str, Sequence[str] | Mapping[str, float]] | None = None,
    subtypes: Mapping[str, Any] | None = None,
    support_weight: float = 2.0,
    anchor_weight: float = 4.0,
    exclusion_weight: float = 2.0,
) -> TypingRegistry:
    """Adapt the calibration registry into broad typing rules.

    Broad marker membership, marker identity/process roles, activity, and SST's
    subtype-only status are read directly from
    :class:`phenocycler.marker_registry.MarkerRegistry`.  With no explicit
    overrides, the versioned :data:`DEFAULT_TYPING_RULES_PATH` supplies the
    curated anchor/support/exclusion and subtype meanings.  The two files have
    independent fingerprints because measurement/calibration policy and
    biological interpretation are separate provenance layers.

    ``anchors``/``exclusions``/``subtypes`` retain a compact explicit-override
    path for tests and pilots. Supplying any of them bypasses the default file;
    they cannot be mixed with ``typing_rules``.
    """
    marker_specs = getattr(marker_registry, "markers", None)
    if marker_specs is None:
        raise TypeError("marker_registry must expose a 'markers' sequence")
    broad: dict[str, dict[str, Any]] = {}
    roles: dict[str, str] = {}

    for marker in marker_specs:
        name = str(getattr(marker, "name"))
        kind = str(getattr(marker, "kind", "type"))
        uses = tuple(str(use) for use in getattr(marker, "uses", ()))
        subtype_only = bool(getattr(marker, "subtype_only", False)) or name in SUBTYPE_ONLY_MARKERS
        if kind == "state":
            roles[name] = PROCESS_ROLE
        elif subtype_only:
            roles[name] = SUBTYPE_ONLY_ROLE
        else:
            roles[name] = IDENTITY_ROLE
        if (
            kind != "type"
            or "broad" not in uses
            or not bool(getattr(marker, "active", False))
        ):
            continue
        compartment = getattr(marker, "broad_compartment", None)
        if compartment is None:
            raise ValueError(f"active broad marker {name!r} has no broad_compartment")
        broad.setdefault(str(compartment), {"markers": []})["markers"].append(name)

    if not broad:
        raise ValueError("marker registry contains no active broad type markers")

    explicit_override = anchors is not None or exclusions is not None or subtypes is not None
    if typing_rules is not None and explicit_override:
        raise ValueError(
            "typing_rules cannot be combined with anchors/exclusions/subtypes overrides"
        )
    if not explicit_override:
        if typing_rules is None:
            overlay = load_typing_rules()
        elif isinstance(typing_rules, TypingRulesOverlay):
            overlay = typing_rules
        elif isinstance(typing_rules, Mapping):
            overlay = typing_rules_from_dict(typing_rules)
        elif isinstance(typing_rules, (str, Path)):
            overlay = load_typing_rules(typing_rules)
        else:
            raise TypeError(
                "typing_rules must be a TypingRulesOverlay, mapping, path, or None"
            )
        validate_typing_rules(overlay, marker_registry)
        rules: dict[str, dict[str, Any]] = {}
        for class_name, rule in overlay.broad.items():
            expanded = dict(rule)
            # Required means donor-marker model validity, not biological
            # indispensability of one protein. Every active model in the
            # registered broad family must be valid-and-negative before Other.
            expanded["required"] = list(broad[class_name]["markers"])
            rules[class_name] = expanded
        return TypingRegistry.from_mapping(
            {
                "marker_roles": roles,
                "broad": rules,
                "subtypes": dict(overlay.subtypes),
                "rules_version": overlay.rules_version,
                "rules_fingerprint": overlay.fingerprint,
                "marker_registry_fingerprint": str(
                    getattr(marker_registry, "fingerprint", "")
                ),
            }
        )

    anchors = dict(anchors or {})
    exclusions = dict(exclusions or {})
    rules: dict[str, dict[str, Any]] = {}
    for class_name, values in broad.items():
        markers = tuple(values["markers"])
        raw_anchors = anchors.get(class_name, ())
        anchor_pairs = _weighted_markers(raw_anchors, default=anchor_weight)
        anchor_names = {marker for marker, _ in anchor_pairs}
        unknown_anchors = sorted(anchor_names - set(markers))
        if unknown_anchors:
            raise ValueError(
                f"{class_name}: anchors are not active broad markers for the class: "
                f"{unknown_anchors}"
            )
        support = {marker: support_weight for marker in markers if marker not in anchor_names}
        raw_exclusions = exclusions.get(class_name, ())
        exclusion_pairs = _weighted_markers(raw_exclusions, default=exclusion_weight)
        rules[class_name] = {
            "anchors": dict(anchor_pairs),
            "support": support,
            "exclusions": dict(exclusion_pairs),
            # Strict Other semantics: every active broad model in the versioned
            # calibration registry must be available and negative.
            "required": list(markers),
        }

    return TypingRegistry.from_mapping(
        {
            "marker_roles": roles,
            "broad": rules,
            "subtypes": dict(subtypes or {}),
        }
    )


def _model_markers(rules: Sequence[ClassRule]) -> tuple[str, ...]:
    return tuple(sorted({marker for rule in rules for marker in rule.model_markers}))


def _evidence_matrices(
    prepared: pd.DataFrame,
    markers: Sequence[str],
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    """Return keys plus aligned support/state/model-valid matrices."""
    markers = tuple(markers)
    keys = prepared.loc[:, list(KEY_COLUMNS)].drop_duplicates(ignore_index=True)
    key_index = pd.MultiIndex.from_frame(keys)

    def pivot(column: str, fill: Any) -> np.ndarray:
        wide = prepared.pivot(index=list(KEY_COLUMNS), columns=MARKER_COLUMN, values=column)
        wide = wide.reindex(index=key_index, columns=markers)
        if fill == "":
            return wide.astype("string").fillna("").to_numpy(dtype=object)
        if fill is False:
            return wide.astype("boolean").fillna(False).to_numpy(dtype=bool)
        return wide.to_numpy()

    support_column = (
        "typing_support" if "typing_support" in prepared else "evidence_probability"
    )
    support = pivot(support_column, None).astype(float)
    state = pivot("state", "")
    valid = pivot("model_valid", False)
    return keys, support, state, valid


def wide_evidence_matrices(
    evidence: pd.DataFrame,
    markers: Sequence[str],
    *,
    tail_probability_is_pvalue: bool = True,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    """Read the production one-row-per-cell calibration schema directly.

    Each marker may provide ``<marker>__log2_threshold_ratio`` (preferred),
    ``__state``, explicit numeric-only ``__evidence_probability``, and
    ``__expression_probability``/``__empirical_tail_probability`` for audit.
    This is the wide equivalent of
    :func:`prepare_cell_marker_evidence`, but it avoids expanding a
    million-cell donor into tens of millions of temporary rows.
    """

    if not isinstance(evidence, pd.DataFrame):
        raise TypeError("evidence must be a pandas DataFrame")
    d = evidence.copy()
    for alias, canonical in {"donor": DONOR_COLUMN, "object": OBJECT_COLUMN}.items():
        if canonical not in d and alias in d:
            d = d.rename(columns={alias: canonical})
    missing = [column for column in KEY_COLUMNS if column not in d]
    if missing:
        raise KeyError(f"wide cell evidence is missing key columns {missing}")
    d[DONOR_COLUMN] = d[DONOR_COLUMN].astype(str)
    d[OBJECT_COLUMN] = d[OBJECT_COLUMN].astype(str)
    if d.duplicated(list(KEY_COLUMNS)).any():
        raise ValueError("wide cell evidence has duplicate donor/object rows")

    keys = d.loc[:, list(KEY_COLUMNS)].reset_index(drop=True)
    n = len(d)
    marker_tuple = tuple(markers)
    support = np.full((n, len(marker_tuple)), np.nan, dtype=float)
    state = np.full((n, len(marker_tuple)), "", dtype=object)
    valid = np.zeros((n, len(marker_tuple)), dtype=bool)

    for index, marker in enumerate(marker_tuple):
        prefix = f"{marker}__"
        state_column = prefix + "state"
        marker_state = (
            _canonical_states(d[state_column])
            if state_column in d
            else pd.Series("", index=d.index, dtype=object)
        )
        contract = _state_safe_marker_support(
            d,
            marker_state,
            ratio_columns=(prefix + "log2_threshold_ratio",),
            score_columns=(prefix + "evidence_probability",),
            tail_columns=(prefix + "empirical_tail_probability",),
            audit_only_columns=(prefix + "expression_probability",),
            explicit_valid_column=prefix + "model_valid",
            tail_probability_is_pvalue=tail_probability_is_pvalue,
            context=marker,
        )
        support[:, index] = contract["typing_support"].to_numpy(dtype=float)
        state[:, index] = marker_state.to_numpy(dtype=object)
        valid[:, index] = contract["model_valid"].to_numpy(dtype=bool)
    return keys, support, state, valid


def evidence_feature_matrix(
    evidence: pd.DataFrame,
    markers: Sequence[str],
    *,
    tail_probability_is_pvalue: bool = True,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    """Return parsed keys and a finite state-safe feature matrix.

    ``features`` is ready for repeated matrix multiplication by fitted or
    bootstrap coefficient matrices: missing, invalid, indeterminate, and
    unavailable entries are all zero. ``state`` and ``valid`` retain the
    categorical decision and per-marker pass mask needed for authoritative
    anchors and strict ``Other``/``Unavailable`` decisions. Long and production
    wide evidence use the same transformer and preserve the requested marker
    order.
    """

    if not isinstance(evidence, pd.DataFrame):
        raise TypeError("evidence must be a pandas DataFrame")
    marker_tuple = tuple(str(marker) for marker in markers)
    if len(set(marker_tuple)) != len(marker_tuple) or any(not marker for marker in marker_tuple):
        raise ValueError("markers must be unique non-empty names")
    if MARKER_COLUMN in evidence:
        prepared = prepare_cell_marker_evidence(
            evidence,
            tail_probability_is_pvalue=tail_probability_is_pvalue,
        )
        keys, support, state, valid = _evidence_matrices(prepared, marker_tuple)
    else:
        keys, support, state, valid = wide_evidence_matrices(
            evidence,
            marker_tuple,
            tail_probability_is_pvalue=tail_probability_is_pvalue,
        )
    features = np.where(valid & np.isfinite(support), support, 0.0)
    if not np.isfinite(features).all() or (features < 0.0).any() or (features > 1.0).any():
        raise AssertionError("state-safe typing features must be finite and lie in [0, 1]")
    return keys, features, state, valid


def balanced_seed_weights(
    seed_labels: pd.DataFrame,
    *,
    class_col: str = "label",
    donor_col: str = DONOR_COLUMN,
) -> pd.Series:
    """Give every class equal mass and every donor within a class equal mass.

    Within one donor/class stratum, its seed cells share that mass equally.
    Returned weights sum to one.  This prevents both common classes and donors
    with many annotated cells from setting the classifier scale.
    """
    if class_col not in seed_labels or donor_col not in seed_labels:
        raise KeyError(f"seed labels require {class_col!r} and {donor_col!r}")
    d = seed_labels[[class_col, donor_col]].copy()
    d[class_col] = d[class_col].astype(str)
    d[donor_col] = d[donor_col].astype(str)
    if d.empty:
        raise ValueError("at least one seed label is required")
    n_classes = d[class_col].nunique()
    donor_counts = d.groupby(class_col, observed=True)[donor_col].transform("nunique")
    cell_counts = d.groupby([class_col, donor_col], observed=True)[donor_col].transform("size")
    weights = 1.0 / (n_classes * donor_counts.astype(float) * cell_counts.astype(float))
    weights.index = seed_labels.index
    # Floating error only; normalization makes the public invariant exact enough
    # for downstream audit tables.
    return weights / weights.sum()


def _finite_softmax(logits: np.ndarray, *, context: str) -> np.ndarray:
    """Compute softmax probabilities and fail closed on numerical corruption."""

    values = np.asarray(logits, dtype=float)
    if values.ndim != 2 or values.shape[1] < 1:
        raise ValueError(f"{context} logits must be a non-empty 2D class matrix")
    if not np.isfinite(values).all():
        raise FloatingPointError(f"{context} produced non-finite logits")
    log_probability = values - logsumexp(values, axis=1, keepdims=True)
    probability = np.exp(log_probability)
    if (
        not np.isfinite(probability).all()
        or (probability < 0).any()
        or (probability > 1).any()
        or not np.allclose(probability.sum(axis=1), 1.0, atol=1e-12)
    ):
        raise FloatingPointError(f"{context} produced invalid posterior probabilities")
    return probability


@dataclass(frozen=True)
class AdditiveMultinomialModel:
    """No-intercept multinomial model with registry-constrained coefficients."""

    classes: tuple[str, ...]
    markers: tuple[str, ...]
    coefficients: np.ndarray
    rules: tuple[ClassRule, ...]
    level: str = "broad"
    parent: str | None = None
    converged: bool = True
    message: str = "registry weights"
    objective: float = np.nan
    n_seeds: int = 0
    tail_probability_is_pvalue: bool = True

    def __post_init__(self) -> None:
        expected = (len(self.classes), len(self.markers))
        if np.asarray(self.coefficients).shape != expected:
            raise ValueError(
                f"coefficient shape {np.asarray(self.coefficients).shape} does not match {expected}"
            )
        if tuple(rule.name for rule in self.rules) != self.classes:
            raise ValueError("model classes and rules must have the same order")

    @property
    def class_priors(self) -> dict[str, float]:
        prior = 1.0 / len(self.classes)
        return {name: prior for name in self.classes}

    def coefficient_table(self) -> pd.DataFrame:
        return pd.DataFrame(self.coefficients, index=self.classes, columns=self.markers)

    def predict_proba(self, evidence: pd.DataFrame) -> pd.DataFrame:
        keys, x, _state, _valid = evidence_feature_matrix(
            evidence,
            self.markers,
            tail_probability_is_pvalue=self.tail_probability_is_pvalue,
        )
        logits = x @ np.asarray(self.coefficients, dtype=float).T
        probs = _finite_softmax(logits, context=f"{self.level} model")
        return pd.concat(
            [keys.reset_index(drop=True), pd.DataFrame(probs, columns=self.classes)],
            axis=1,
        )


def classifier_from_registry(
    registry: Mapping[str, Any] | TypingRegistry,
    *,
    level: str = "broad",
    parent: str | None = None,
    max_abs_weight: float = 8.0,
    tail_probability_is_pvalue: bool = True,
) -> AdditiveMultinomialModel:
    """Build the transparent, unfitted additive classifier encoded by the registry."""
    reg = _as_registry(registry)
    rules = reg.rules_for(level=level, parent=parent)
    if not rules:
        raise ValueError(f"no {level} rules available" + (f" for parent {parent!r}" if parent else ""))
    markers = _model_markers(rules)
    coefficients = np.zeros((len(rules), len(markers)), dtype=float)
    for ci, rule in enumerate(rules):
        for mi, marker in enumerate(markers):
            coefficients[ci, mi] = np.clip(
                rule.initial_coefficient(marker),
                -max_abs_weight,
                max_abs_weight,
            )
    return AdditiveMultinomialModel(
        classes=tuple(rule.name for rule in rules),
        markers=markers,
        coefficients=coefficients,
        rules=rules,
        level=level,
        parent=parent,
        tail_probability_is_pvalue=tail_probability_is_pvalue,
    )


def _normalise_seed_labels(
    seed_labels: pd.DataFrame,
    *,
    label_col: str,
) -> pd.DataFrame:
    if not isinstance(seed_labels, pd.DataFrame):
        raise TypeError("seed_labels must be a pandas DataFrame")
    d = seed_labels.copy()
    if DONOR_COLUMN not in d and "donor" in d:
        d = d.rename(columns={"donor": DONOR_COLUMN})
    if OBJECT_COLUMN not in d and "object" in d:
        d = d.rename(columns={"object": OBJECT_COLUMN})
    missing = [c for c in (*KEY_COLUMNS, label_col) if c not in d]
    if missing:
        raise KeyError(f"seed labels are missing required columns {missing}")
    d[DONOR_COLUMN] = d[DONOR_COLUMN].astype(str)
    d[OBJECT_COLUMN] = d[OBJECT_COLUMN].astype(str)
    d[label_col] = d[label_col].astype(str)
    conflicting = (
        d.groupby(list(KEY_COLUMNS), observed=True)[label_col].nunique().gt(1)
    )
    if conflicting.any():
        bad = conflicting[conflicting].index.tolist()[:5]
        raise ValueError(f"seed cells have conflicting class labels: {bad}")
    return d


def fit_constrained_classifier(
    evidence: pd.DataFrame,
    registry: Mapping[str, Any] | TypingRegistry,
    seed_labels: pd.DataFrame,
    *,
    level: str = "broad",
    parent: str | None = None,
    label_col: str = "label",
    l2: float = 0.02,
    max_abs_weight: float = 8.0,
    maxiter: int = 500,
    balance_group_col: str | None = None,
    tail_probability_is_pvalue: bool = True,
) -> AdditiveMultinomialModel:
    """Fit a sign-constrained additive multinomial model to curated seed cells.

    There is deliberately no intercept: a cell with no positive marker evidence
    receives the uniform prior.  Anchor/support coefficients are bounded to
    ``[0, max_abs_weight]``; exclusions to ``[-max_abs_weight, 0]``; all
    unlisted marker/class coefficients are fixed at zero.
    """
    if l2 < 0 or not np.isfinite(l2):
        raise ValueError("l2 must be finite and >= 0")
    if not np.isfinite(max_abs_weight) or max_abs_weight <= 0:
        raise ValueError("max_abs_weight must be finite and > 0")

    reg = _as_registry(registry)
    initial = classifier_from_registry(
        reg,
        level=level,
        parent=parent,
        max_abs_weight=max_abs_weight,
        tail_probability_is_pvalue=tail_probability_is_pvalue,
    )
    keys, features, _state, _valid = evidence_feature_matrix(
        evidence,
        initial.markers,
        tail_probability_is_pvalue=tail_probability_is_pvalue,
    )
    feature_frame = pd.concat(
        [
            keys.reset_index(drop=True),
            pd.DataFrame(features, columns=initial.markers),
        ],
        axis=1,
    )
    seeds = _normalise_seed_labels(seed_labels, label_col=label_col)
    unknown = sorted(set(seeds[label_col]) - set(initial.classes))
    if unknown:
        raise ValueError(f"seed labels are not registry classes: {unknown}")
    merged = seeds.merge(feature_frame, on=list(KEY_COLUMNS), how="left", validate="m:1", indicator=True)
    absent = merged["_merge"] != "both"
    if absent.any():
        examples = merged.loc[absent, list(KEY_COLUMNS)].head(5).to_dict("records")
        raise ValueError(f"seed cells have no marker evidence: {examples}")
    merged = merged.drop(columns="_merge")

    balance_col = balance_group_col or DONOR_COLUMN
    if balance_col not in merged:
        raise KeyError(f"balance group column {balance_col!r} is missing from seed labels")
    sample_weight = balanced_seed_weights(
        merged,
        class_col=label_col,
        donor_col=balance_col,
    ).to_numpy(float)
    x = merged.loc[:, list(initial.markers)].to_numpy(float)
    class_index = {name: i for i, name in enumerate(initial.classes)}
    y = merged[label_col].map(class_index).to_numpy(int)

    free: list[tuple[int, int]] = []
    bounds: list[tuple[float, float]] = []
    x0: list[float] = []
    for ci, rule in enumerate(initial.rules):
        positive = set(rule.anchor_markers + rule.support_markers)
        negative = set(rule.exclusion_markers)
        for mi, marker in enumerate(initial.markers):
            if marker in positive:
                free.append((ci, mi))
                bounds.append((0.0, max_abs_weight))
                x0.append(max(0.0, initial.coefficients[ci, mi]))
            elif marker in negative:
                free.append((ci, mi))
                bounds.append((-max_abs_weight, 0.0))
                x0.append(min(0.0, initial.coefficients[ci, mi]))

    base = np.zeros_like(initial.coefficients, dtype=float)

    def unpack(theta: np.ndarray) -> np.ndarray:
        weights = base.copy()
        for value, (ci, mi) in zip(theta, free):
            weights[ci, mi] = value
        return weights

    def objective(theta: np.ndarray) -> tuple[float, np.ndarray]:
        weights = unpack(theta)
        logits = x @ weights.T
        log_prob = logits - logsumexp(logits, axis=1, keepdims=True)
        loss = -float(np.sum(sample_weight * log_prob[np.arange(len(y)), y]))
        prob = np.exp(log_prob)
        residual = prob
        residual[np.arange(len(y)), y] -= 1.0
        grad_full = (sample_weight[:, None] * residual).T @ x
        gradient = np.array([grad_full[ci, mi] for ci, mi in free], dtype=float)
        if l2:
            loss += 0.5 * l2 * float(theta @ theta)
            gradient += l2 * theta
        return loss, gradient

    if free:
        result = minimize(
            objective,
            np.asarray(x0, dtype=float),
            method="L-BFGS-B",
            jac=True,
            bounds=bounds,
            options={"maxiter": int(maxiter)},
        )
        coefficients = unpack(np.asarray(result.x, dtype=float))
        converged = bool(result.success)
        message = str(result.message)
        value = float(result.fun)
    else:
        coefficients = base
        converged = True
        message = "no free coefficients"
        value = float(objective(np.empty(0))[0])

    return AdditiveMultinomialModel(
        classes=initial.classes,
        markers=initial.markers,
        coefficients=coefficients,
        rules=initial.rules,
        level=level,
        parent=parent,
        converged=converged,
        message=message,
        objective=value,
        n_seeds=len(merged),
        tail_probability_is_pvalue=tail_probability_is_pvalue,
    )


@dataclass(frozen=True)
class AssignmentThresholds:
    """Frozen acceptance thresholds for inferred identity calls."""

    inferred_probability: float = 0.80
    inferred_margin: float = 0.25
    inferred_stability: float = 0.90
    anchor_probability: float = 0.95
    negative_probability: float = 0.20

    def __post_init__(self) -> None:
        values = {
            "inferred_probability": self.inferred_probability,
            "inferred_margin": self.inferred_margin,
            "inferred_stability": self.inferred_stability,
            "anchor_probability": self.anchor_probability,
            "negative_probability": self.negative_probability,
        }
        bad = {name: value for name, value in values.items() if not 0 <= value <= 1}
        if bad:
            raise ValueError(f"assignment thresholds must lie in [0, 1]: {bad}")


def _align_stability(keys: pd.DataFrame, stability: Any, *, prefix: str) -> np.ndarray:
    if stability is None:
        return np.full(len(keys), np.nan)
    if np.isscalar(stability):
        value = float(stability)
        if not 0 <= value <= 1:
            raise ValueError("stability must lie in [0, 1]")
        return np.full(len(keys), value)
    if isinstance(stability, pd.Series):
        if isinstance(stability.index, pd.MultiIndex):
            lookup = stability.rename("stability").reset_index()
        elif len(stability) == len(keys):
            values = stability.to_numpy(float)
            if np.any((values < 0) | (values > 1)):
                raise ValueError("stability must lie in [0, 1]")
            return values
        else:
            raise ValueError("stability Series must align by a donor/object MultiIndex or row count")
    elif isinstance(stability, pd.DataFrame):
        lookup = stability.copy()
    else:
        raise TypeError("stability must be None, a scalar, Series, or DataFrame")
    if DONOR_COLUMN not in lookup and "donor" in lookup:
        lookup = lookup.rename(columns={"donor": DONOR_COLUMN})
    if OBJECT_COLUMN not in lookup and "object" in lookup:
        lookup = lookup.rename(columns={"object": OBJECT_COLUMN})
    value_columns = [
        c for c in ("stability", f"{prefix}_stability") if c in lookup
    ]
    if not value_columns or any(c not in lookup for c in KEY_COLUMNS):
        raise KeyError(
            f"stability table requires {list(KEY_COLUMNS)} and 'stability' or "
            f"{prefix + '_stability'!r}"
        )
    value_col = value_columns[0]
    lookup[DONOR_COLUMN] = lookup[DONOR_COLUMN].astype(str)
    lookup[OBJECT_COLUMN] = lookup[OBJECT_COLUMN].astype(str)
    if lookup.duplicated(list(KEY_COLUMNS)).any():
        raise ValueError("stability table has duplicate donor/object rows")
    aligned = keys.merge(
        lookup[[*KEY_COLUMNS, value_col]],
        on=list(KEY_COLUMNS),
        how="left",
        validate="1:1",
    )[value_col].to_numpy(float)
    if np.any((aligned[np.isfinite(aligned)] < 0) | (aligned[np.isfinite(aligned)] > 1)):
        raise ValueError("stability must lie in [0, 1]")
    return aligned


def _positive_marker(
    probability: np.ndarray,
    state: np.ndarray,
    valid: np.ndarray,
    marker_index: int,
    threshold: float,
) -> np.ndarray:
    # A categorical calibration call is authoritative when present.  In
    # particular, a valid ``negative`` call must not be promoted to an anchor
    # merely because its continuous ranking probability is high but still
    # short of that marker's stricter familywise alpha.  Numeric-only pilot
    # evidence retains the explicit threshold fallback.
    return valid[:, marker_index] & (
        (state[:, marker_index] == "positive")
        | (
            (state[:, marker_index] == "")
            & (probability[:, marker_index] >= threshold)
        )
    )


def _assign_level(
    evidence: pd.DataFrame,
    model: AdditiveMultinomialModel,
    *,
    stability: Any,
    thresholds: AssignmentThresholds,
    prefix: str,
) -> pd.DataFrame:
    if MARKER_COLUMN in evidence:
        prepared = prepare_cell_marker_evidence(
            evidence,
            tail_probability_is_pvalue=model.tail_probability_is_pvalue,
        )
        keys, probability, state, valid = _evidence_matrices(
            prepared, model.markers
        )
    else:
        keys, probability, state, valid = wide_evidence_matrices(
            evidence,
            model.markers,
            tail_probability_is_pvalue=model.tail_probability_is_pvalue,
        )
    marker_index = {marker: i for i, marker in enumerate(model.markers)}
    x = np.where(valid & np.isfinite(probability), probability, 0.0)
    valid_positive_support = np.zeros((len(keys), len(model.classes)), dtype=bool)
    for ci, rule in enumerate(model.rules):
        for marker in dict.fromkeys(rule.anchor_markers + rule.support_markers):
            mi = marker_index[marker]
            valid_positive_support[:, ci] |= valid[:, mi] & (
                (state[:, mi] == "positive") | (x[:, mi] > 0.0)
            )
    logits = x @ np.asarray(model.coefficients, dtype=float).T
    posterior = _finite_softmax(logits, context=f"{prefix} assignment model")
    n = len(keys)

    evidence_exists = np.isfinite(probability).any(axis=1)
    order = np.argsort(-logits, axis=1, kind="stable")
    best_index = order[:, 0]
    best_probability = posterior[np.arange(n), best_index]
    if len(model.classes) > 1:
        second_probability = posterior[np.arange(n), order[:, 1]]
        winner_logit_gap = (
            logits[np.arange(n), best_index]
            - logits[np.arange(n), order[:, 1]]
        )
        inference_unique_winner_pass = winner_logit_gap > 1e-12
    else:
        second_probability = np.zeros(n, dtype=float)
        inference_unique_winner_pass = np.ones(n, dtype=bool)
    margin = best_probability - second_probability
    best_candidate = np.asarray(model.classes, dtype=object)[best_index]
    best_candidate = np.where(evidence_exists, best_candidate, "")
    best_probability = np.where(evidence_exists, best_probability, np.nan)
    margin = np.where(evidence_exists, margin, np.nan)
    aligned_stability = _align_stability(keys, stability, prefix=prefix)
    inference_probability_pass = (
        np.isfinite(best_probability)
        & (best_probability >= thresholds.inferred_probability)
    )
    inference_margin_pass = np.isfinite(margin) & (
        margin >= thresholds.inferred_margin
    )
    inference_stability_pass = np.isfinite(aligned_stability) & (
        aligned_stability >= thresholds.inferred_stability
    )
    inference_valid_evidence_pass = (
        evidence_exists
        & valid_positive_support[np.arange(n), best_index]
    )

    anchor_hits: list[list[str]] = [[] for _ in range(n)]
    anchor_classes: list[list[str]] = [[] for _ in range(n)]
    exclusion_conflict = np.zeros((n, len(model.classes)), dtype=bool)
    for ci, rule in enumerate(model.rules):
        class_anchor = np.zeros(n, dtype=bool)
        for marker in rule.anchor_markers:
            mi = marker_index[marker]
            hit = _positive_marker(
                probability,
                state,
                valid,
                mi,
                thresholds.anchor_probability,
            )
            class_anchor |= hit
            for row in np.flatnonzero(hit):
                anchor_hits[row].append(marker)
        for row in np.flatnonzero(class_anchor):
            anchor_classes[row].append(rule.name)
        for marker in rule.exclusion_markers:
            mi = marker_index[marker]
            exclusion_conflict[:, ci] |= _positive_marker(
                probability,
                state,
                valid,
                mi,
                thresholds.anchor_probability,
            )

    all_required_valid = np.ones(n, dtype=bool)
    all_required_negative = np.ones(n, dtype=bool)
    for rule in model.rules:
        for marker in rule.required:
            mi = marker_index[marker]
            marker_valid = valid[:, mi] & np.isfinite(probability[:, mi])
            marker_negative = marker_valid & (
                (state[:, mi] == "negative")
                | (
                    (state[:, mi] != "positive")
                    & (probability[:, mi] <= thresholds.negative_probability)
                )
            )
            all_required_valid &= marker_valid
            all_required_negative &= marker_negative
    all_required_negative &= all_required_valid

    positive_evidence = np.zeros(n, dtype=bool)
    for mi in range(len(model.markers)):
        positive_evidence |= valid[:, mi] & (
            (state[:, mi] == "positive")
            | (
                (state[:, mi] == "")
                & (probability[:, mi] > thresholds.negative_probability)
            )
        )

    label = np.full(n, UNAVAILABLE_LABEL, dtype=object)
    status = np.full(n, UNAVAILABLE, dtype=object)
    reason = np.full(n, "no_identity_evidence", dtype=object)

    for row in range(n):
        if not evidence_exists[row]:
            continue
        anchors = list(dict.fromkeys(anchor_classes[row]))
        if len(anchors) == 1:
            anchor_class = anchors[0]
            ci = model.classes.index(anchor_class)
            if not exclusion_conflict[row, ci]:
                label[row] = anchor_class
                best_candidate[row] = anchor_class
                best_probability[row] = posterior[row, ci]
                if len(model.classes) > 1:
                    margin[row] = posterior[row, ci] - np.max(
                        np.delete(posterior[row], ci)
                    )
                else:
                    margin[row] = posterior[row, ci]
                status[row] = AUTHORITATIVE_ANCHOR
                reason[row] = "single_authoritative_anchor"
                continue
            label[row] = AMBIGUOUS_LABEL
            status[row] = AMBIGUOUS
            reason[row] = "anchor_exclusion_conflict"
            continue
        if len(anchors) > 1:
            label[row] = AMBIGUOUS_LABEL
            status[row] = AMBIGUOUS
            reason[row] = "conflicting_authoritative_anchors"
            continue
        if all_required_negative[row]:
            label[row] = OTHER_LABEL
            status[row] = OTHER
            reason[row] = "all_required_models_valid_and_negative"
            continue
        if (
            inference_unique_winner_pass[row]
            and inference_probability_pass[row]
            and inference_margin_pass[row]
            and inference_stability_pass[row]
            and inference_valid_evidence_pass[row]
        ):
            label[row] = best_candidate[row]
            status[row] = INFERRED
            reason[row] = "posterior_margin_stability_pass"
            continue
        if not positive_evidence[row] and not all_required_valid[row]:
            label[row] = UNAVAILABLE_LABEL
            status[row] = UNAVAILABLE
            reason[row] = "required_models_unavailable"
            continue
        label[row] = AMBIGUOUS_LABEL
        status[row] = AMBIGUOUS
        reason[row] = (
            "nonunique_model_winner"
            if not inference_unique_winner_pass[row]
            else (
                "valid_positive_support_required"
                if not inference_valid_evidence_pass[row]
                else "posterior_margin_or_stability_failed"
            )
        )

    if np.any(evidence_exists & (best_candidate == "")):
        raise AssertionError("a best candidate is required wherever identity evidence exists")

    # Recompute the public pass diagnostics after authoritative-anchor handling,
    # which can intentionally replace the posterior's original best candidate.
    class_index = {class_name: index for index, class_name in enumerate(model.classes)}
    final_best_index = best_index.copy()
    for row in np.flatnonzero(evidence_exists):
        final_best_index[row] = class_index[str(best_candidate[row])]
    probability_pass = np.isfinite(best_probability) & (
        best_probability >= thresholds.inferred_probability
    )
    margin_pass = np.isfinite(margin) & (margin >= thresholds.inferred_margin)
    stability_pass = np.isfinite(aligned_stability) & (
        aligned_stability >= thresholds.inferred_stability
    )
    valid_evidence_pass = (
        evidence_exists
        & valid_positive_support[np.arange(n), final_best_index]
    )
    threshold_eligible = np.asarray(
        [
            evidence_exists[row]
            and not list(dict.fromkeys(anchor_classes[row]))
            and not all_required_negative[row]
            and not (
                not positive_evidence[row]
                and not all_required_valid[row]
            )
            for row in range(n)
        ],
        dtype=bool,
    )
    inference_failure_reasons: list[str] = []
    for row in range(n):
        failures: list[str] = []
        if not evidence_exists[row]:
            failures.append("no_identity_evidence")
        else:
            if not inference_unique_winner_pass[row]:
                failures.append("best_candidate_not_unique")
            if not valid_evidence_pass[row]:
                failures.append("no_valid_positive_support")
            if not probability_pass[row]:
                failures.append("probability_below_threshold")
            if not margin_pass[row]:
                failures.append("margin_below_threshold")
            if not np.isfinite(aligned_stability[row]):
                failures.append("stability_unavailable")
            elif not stability_pass[row]:
                failures.append("stability_below_threshold")
        inference_failure_reasons.append("|".join(failures))

    out = keys.copy()
    out[f"{prefix}_label"] = label
    out[f"{prefix}_assignment_status"] = status
    out[f"{prefix}_best_candidate"] = best_candidate
    out[f"{prefix}_best_probability"] = best_probability
    out[f"{prefix}_margin"] = margin
    out[f"{prefix}_stability"] = aligned_stability
    out[f"{prefix}_unique_winner_pass"] = inference_unique_winner_pass
    out[f"{prefix}_probability_pass"] = probability_pass
    out[f"{prefix}_margin_pass"] = margin_pass
    out[f"{prefix}_stability_pass"] = stability_pass
    out[f"{prefix}_valid_evidence_pass"] = valid_evidence_pass
    out[f"{prefix}_threshold_eligible"] = threshold_eligible
    out[f"{prefix}_inference_failure_reasons"] = inference_failure_reasons
    out[f"{prefix}_reason"] = reason
    out[f"{prefix}_all_required_valid"] = all_required_valid
    out[f"{prefix}_all_required_negative"] = all_required_negative
    out[f"{prefix}_authoritative_markers"] = [
        "|".join(dict.fromkeys(markers)) for markers in anchor_hits
    ]
    out[f"{prefix}_anchor_classes"] = [
        "|".join(dict.fromkeys(classes)) for classes in anchor_classes
    ]
    for ci, class_name in enumerate(model.classes):
        out[f"{prefix}_probability__{class_name}"] = posterior[:, ci]
    rankings = []
    for row in range(n):
        if not evidence_exists[row]:
            rankings.append("[]")
            continue
        ranked = sorted(
            (
                (class_name, float(posterior[row, ci]))
                for ci, class_name in enumerate(model.classes)
            ),
            key=lambda item: (-item[1], model.classes.index(item[0])),
        )
        rankings.append(json.dumps(ranked, separators=(",", ":")))
    out[f"ranked_{prefix}_probabilities"] = rankings
    if prefix == "broad":
        # Stable export-facing contract; detailed broad_* fields above remain
        # available for audits and downstream subtype decisions.
        out["broad_type"] = out["broad_label"]
        out["best_broad_type"] = out["broad_best_candidate"]
        out["assignment_status"] = out["broad_assignment_status"]
        out["confidence"] = out["broad_best_probability"]
        out["ranked_probabilities"] = out["ranked_broad_probabilities"]
    return out


def assign_broad_types(
    evidence: pd.DataFrame,
    registry: Mapping[str, Any] | TypingRegistry,
    *,
    model: AdditiveMultinomialModel | None = None,
    stability: Any = None,
    thresholds: AssignmentThresholds | None = None,
) -> pd.DataFrame:
    """Assign simultaneous broad identity under the authoritative output contract.

    Non-anchor inference requires probability >= 0.80, margin >= 0.25, and
    grouped-donor stability >= 0.90 by default.  ``Other`` is emitted only when
    every required broad marker model is available and negative.
    """
    reg = _as_registry(registry)
    model = model or classifier_from_registry(reg, level="broad")
    if model.level != "broad":
        raise ValueError("assign_broad_types requires a broad model")
    return _assign_level(
        evidence,
        model,
        stability=stability,
        thresholds=thresholds or AssignmentThresholds(),
        prefix="broad",
    )


def assign_subtypes(
    evidence: pd.DataFrame,
    broad_assignments: pd.DataFrame,
    registry: Mapping[str, Any] | TypingRegistry,
    *,
    models: Mapping[str, AdditiveMultinomialModel] | None = None,
    stability: Mapping[str, Any] | None = None,
    thresholds: AssignmentThresholds | None = None,
) -> pd.DataFrame:
    """Apply subtype models only inside an accepted parent broad class.

    A parent with configured subtypes but no accepted subtype is reported as
    ``<Parent>-unclassified``.  A terminal parent with no subtype rules retains
    its broad label.
    """
    reg = _as_registry(registry)
    thresholds = thresholds or AssignmentThresholds()
    models = dict(models or {})
    stability = dict(stability or {})
    required_broad = [*KEY_COLUMNS, "broad_label", "broad_assignment_status"]
    missing = [column for column in required_broad if column not in broad_assignments]
    if missing:
        raise KeyError(f"broad assignments are missing required columns {missing}")
    if broad_assignments.duplicated(list(KEY_COLUMNS)).any():
        raise ValueError("broad assignments have duplicate donor/object rows")

    out = broad_assignments.copy()
    out[DONOR_COLUMN] = out[DONOR_COLUMN].astype(str)
    out[OBJECT_COLUMN] = out[OBJECT_COLUMN].astype(str)
    out["subtype_label"] = ""
    out["subtype_assignment_status"] = UNAVAILABLE
    out["subtype_best_candidate"] = ""
    out["subtype_best_probability"] = np.nan
    out["subtype_margin"] = np.nan
    out["subtype_stability"] = np.nan
    out["subtype_unique_winner_pass"] = False
    out["subtype_probability_pass"] = False
    out["subtype_margin_pass"] = False
    out["subtype_stability_pass"] = False
    out["subtype_valid_evidence_pass"] = False
    out["subtype_threshold_eligible"] = False
    out["subtype_inference_failure_reasons"] = "parent_not_classified"
    out["subtype_reason"] = "parent_not_classified"
    out["ranked_subtype_probabilities"] = "[]"
    out["cell_type"] = out["broad_label"]

    if MARKER_COLUMN in evidence:
        evidence_keys = _resolve_key_columns(evidence)
    else:
        evidence_keys = evidence.copy()
        for alias, canonical in {"donor": DONOR_COLUMN, "object": OBJECT_COLUMN}.items():
            if canonical not in evidence_keys and alias in evidence_keys:
                evidence_keys = evidence_keys.rename(columns={alias: canonical})
        missing_evidence_keys = [
            column for column in KEY_COLUMNS if column not in evidence_keys
        ]
        if missing_evidence_keys:
            raise KeyError(
                f"wide cell evidence is missing key columns {missing_evidence_keys}"
            )
        if evidence_keys.duplicated(list(KEY_COLUMNS)).any():
            raise ValueError("wide cell evidence has duplicate donor/object rows")
    evidence_keys[DONOR_COLUMN] = evidence_keys[DONOR_COLUMN].astype(str)
    evidence_keys[OBJECT_COLUMN] = evidence_keys[OBJECT_COLUMN].astype(str)

    accepted_broad = out["broad_assignment_status"].isin(
        {AUTHORITATIVE_ANCHOR, INFERRED}
    )
    for parent in (rule.name for rule in reg.broad_rules):
        parent_mask = accepted_broad & (out["broad_label"] == parent)
        if not parent_mask.any():
            continue
        rules = reg.rules_for(level="subtype", parent=parent)
        if not rules:
            out.loc[parent_mask, "subtype_assignment_status"] = "not_applicable"
            out.loc[parent_mask, "subtype_reason"] = "terminal_parent"
            out.loc[
                parent_mask, "subtype_inference_failure_reasons"
            ] = "not_applicable"
            out.loc[parent_mask, "cell_type"] = parent
            continue

        parent_keys = out.loc[parent_mask, list(KEY_COLUMNS)]
        parent_evidence = evidence_keys.merge(
            parent_keys,
            on=list(KEY_COLUMNS),
            how="inner",
            validate="m:1",
        )
        parent_model = models.get(parent) or classifier_from_registry(
            reg,
            level="subtype",
            parent=parent,
        )
        if parent_model.level != "subtype" or parent_model.parent != parent:
            raise ValueError(f"subtype model for {parent!r} has the wrong level or parent")
        assigned = _assign_level(
            parent_evidence,
            parent_model,
            stability=stability.get(parent),
            thresholds=thresholds,
            prefix="subtype",
        )
        assigned = assigned[
            [
                *KEY_COLUMNS,
                "subtype_label",
                "subtype_assignment_status",
                "subtype_best_candidate",
                "subtype_best_probability",
                "subtype_margin",
                "subtype_stability",
                "subtype_unique_winner_pass",
                "subtype_probability_pass",
                "subtype_margin_pass",
                "subtype_stability_pass",
                "subtype_valid_evidence_pass",
                "subtype_threshold_eligible",
                "subtype_inference_failure_reasons",
                "subtype_reason",
                "ranked_subtype_probabilities",
            ]
        ]
        index = pd.MultiIndex.from_frame(out[list(KEY_COLUMNS)])
        assigned_index = pd.MultiIndex.from_frame(assigned[list(KEY_COLUMNS)])
        positions = index.get_indexer(assigned_index)
        if (positions < 0).any():
            raise AssertionError("subtype assignments escaped their parent cell set")
        for column in assigned.columns.difference(KEY_COLUMNS, sort=False):
            out.iloc[positions, out.columns.get_loc(column)] = assigned[column].to_numpy()

        accepted_subtype = out.loc[parent_mask, "subtype_assignment_status"].isin(
            {AUTHORITATIVE_ANCHOR, INFERRED}
        )
        parent_rows = out.index[parent_mask]
        accepted_rows = parent_rows[accepted_subtype.to_numpy()]
        fallback_rows = parent_rows[~accepted_subtype.to_numpy()]
        out.loc[accepted_rows, "cell_type"] = out.loc[accepted_rows, "subtype_label"]
        out.loc[fallback_rows, "cell_type"] = f"{parent}-unclassified"
        out.loc[fallback_rows, "subtype_label"] = f"{parent}-unclassified"

    # Final export aliases describe the most specific accepted/fallback result,
    # while broad_* and subtype_* retain the two decisions separately.
    out["broad_type"] = out["broad_label"]
    out["specific_type"] = out["cell_type"]
    out["best_broad_type"] = out["broad_best_candidate"]
    out["assignment_status"] = out["broad_assignment_status"]
    out["confidence"] = out["broad_best_probability"]
    out["best_specific_type"] = out["broad_best_candidate"]
    out["ranked_specific_probabilities"] = out[
        "ranked_broad_probabilities"
    ]
    subtype_evaluated = out["subtype_assignment_status"] != "not_applicable"
    accepted_parent = out["broad_assignment_status"].isin(
        {AUTHORITATIVE_ANCHOR, INFERRED}
    )
    use_subtype_status = subtype_evaluated & accepted_parent
    out.loc[use_subtype_status, "assignment_status"] = out.loc[
        use_subtype_status, "subtype_assignment_status"
    ]
    accepted_specific = out["subtype_assignment_status"].isin(
        {AUTHORITATIVE_ANCHOR, INFERRED}
    )
    both_confidence = np.minimum(
        pd.to_numeric(out["broad_best_probability"], errors="coerce"),
        pd.to_numeric(out["subtype_best_probability"], errors="coerce"),
    )
    out.loc[accepted_specific, "confidence"] = both_confidence[accepted_specific]
    has_subtype_ranking = (
        accepted_parent
        & subtype_evaluated
        & out["ranked_subtype_probabilities"].fillna("[]").ne("[]")
    )
    out.loc[has_subtype_ranking, "best_specific_type"] = out.loc[
        has_subtype_ranking, "subtype_best_candidate"
    ]
    out.loc[
        has_subtype_ranking, "ranked_specific_probabilities"
    ] = out.loc[has_subtype_ranking, "ranked_subtype_probabilities"]
    out["specific_assignment_status"] = out["assignment_status"]
    # Backward-compatible generic ranking now means the most specific level
    # actually evaluated; broad alternatives remain explicit in
    # ranked_broad_probabilities.
    out["ranked_probabilities"] = out["ranked_specific_probabilities"]
    return out


def assign_hierarchical_types(
    evidence: pd.DataFrame,
    registry: Mapping[str, Any] | TypingRegistry,
    *,
    broad_model: AdditiveMultinomialModel | None = None,
    subtype_models: Mapping[str, AdditiveMultinomialModel] | None = None,
    broad_stability: Any = None,
    subtype_stability: Mapping[str, Any] | None = None,
    thresholds: AssignmentThresholds | None = None,
) -> pd.DataFrame:
    """Run broad assignment followed by parent-constrained subtype assignment."""
    reg = _as_registry(registry)
    broad = assign_broad_types(
        evidence,
        reg,
        model=broad_model,
        stability=broad_stability,
        thresholds=thresholds,
    )
    return assign_subtypes(
        evidence,
        broad,
        reg,
        models=subtype_models,
        stability=subtype_stability,
        thresholds=thresholds,
    )


def grouped_donor_bootstrap_stability(
    evidence: pd.DataFrame,
    registry: Mapping[str, Any] | TypingRegistry,
    seed_labels: pd.DataFrame,
    *,
    model: AdditiveMultinomialModel | None = None,
    level: str = "broad",
    parent: str | None = None,
    label_col: str = "label",
    n_bootstrap: int = 100,
    random_state: int | np.random.Generator | None = 0,
    l2: float = 0.02,
    max_abs_weight: float = 8.0,
    maxiter: int = 500,
    tail_probability_is_pvalue: bool = True,
) -> pd.DataFrame:
    """Estimate best-candidate stability by resampling whole seed donors.

    Each bootstrap draw is treated as a distinct balancing group, so drawing a
    donor twice has the intended multiplicity while cells within each
    donor/class draw remain balanced.
    """
    if int(n_bootstrap) != n_bootstrap or n_bootstrap < 1:
        raise ValueError("n_bootstrap must be an integer >= 1")
    reg = _as_registry(registry)
    seeds = _normalise_seed_labels(seed_labels, label_col=label_col)
    donors = seeds[DONOR_COLUMN].drop_duplicates().to_numpy(object)
    if len(donors) == 0:
        raise ValueError("at least one seed donor is required")
    if isinstance(random_state, np.random.Generator):
        rng = random_state
    else:
        rng = np.random.default_rng(random_state)
    model = model or fit_constrained_classifier(
        evidence,
        reg,
        seeds,
        level=level,
        parent=parent,
        label_col=label_col,
        l2=l2,
        max_abs_weight=max_abs_weight,
        maxiter=maxiter,
        tail_probability_is_pvalue=tail_probability_is_pvalue,
    )
    base = model.predict_proba(evidence)
    base_values = base.loc[:, list(model.classes)].to_numpy(float)
    base_best = np.asarray(model.classes, dtype=object)[np.argmax(base_values, axis=1)]
    matches = np.zeros(len(base), dtype=int)
    completed = 0

    for _ in range(int(n_bootstrap)):
        draw = rng.choice(donors, size=len(donors), replace=True)
        parts: list[pd.DataFrame] = []
        for draw_index, donor in enumerate(draw):
            part = seeds.loc[seeds[DONOR_COLUMN] == str(donor)].copy()
            part["_bootstrap_group"] = f"{draw_index}:{donor}"
            parts.append(part)
        boot_seeds = pd.concat(parts, ignore_index=True)
        boot_model = fit_constrained_classifier(
            evidence,
            reg,
            boot_seeds,
            level=level,
            parent=parent,
            label_col=label_col,
            l2=l2,
            max_abs_weight=max_abs_weight,
            maxiter=maxiter,
            balance_group_col="_bootstrap_group",
            tail_probability_is_pvalue=tail_probability_is_pvalue,
        )
        boot_prob = boot_model.predict_proba(evidence).loc[:, list(model.classes)].to_numpy(float)
        boot_best = np.asarray(model.classes, dtype=object)[np.argmax(boot_prob, axis=1)]
        matches += boot_best == base_best
        completed += 1

    out = base.loc[:, list(KEY_COLUMNS)].copy()
    out["best_candidate"] = base_best
    out["stability"] = matches / completed
    out["bootstrap_replicates"] = completed
    return out


__all__ = [
    "AMBIGUOUS",
    "AMBIGUOUS_LABEL",
    "AUTHORITATIVE_ANCHOR",
    "AdditiveMultinomialModel",
    "AssignmentThresholds",
    "ClassRule",
    "DEFAULT_TYPING_RULES_PATH",
    "FEATURE_SCHEMA_VERSION",
    "IDENTITY_ROLE",
    "INFERRED",
    "OTHER",
    "OTHER_LABEL",
    "PROCESS_ROLE",
    "SUBTYPE_ONLY_MARKERS",
    "SUBTYPE_ONLY_ROLE",
    "TypingRegistry",
    "TypingRulesOverlay",
    "UNAVAILABLE",
    "UNAVAILABLE_LABEL",
    "assign_broad_types",
    "assign_hierarchical_types",
    "assign_subtypes",
    "balanced_seed_weights",
    "classifier_from_registry",
    "evidence_feature_matrix",
    "fit_constrained_classifier",
    "grouped_donor_bootstrap_stability",
    "load_typing_rules",
    "prepare_cell_marker_evidence",
    "typing_rules_from_dict",
    "typing_registry_from_marker_registry",
    "wide_evidence_matrices",
    "validate_typing_rules",
]
