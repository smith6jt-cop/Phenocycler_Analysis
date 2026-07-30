"""Validated marker registry for donor-local expression calibration.

The JSON beside this module is deliberately the source of truth.  This module
only supplies typed access, validation, and deterministic familywise error
allocation; it does not copy marker lists from ``config`` or
``marker_taxonomy``.

An ``active`` marker has one pre-declared mutually exclusive reference.  A
marker can still be present as ``deferred`` or ``standalone`` so that absence
of a calibrated call is explicit rather than silently interpreted as a
negative result.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping


DEFAULT_REGISTRY_PATH = Path(__file__).with_name("marker_registry.json")
ALLOWED_KINDS = frozenset({"type", "state", "structural", "excluded"})
ALLOWED_CALIBRATION_STATUSES = frozenset({"active", "deferred", "standalone"})
ALLOWED_BROAD_COMPARTMENTS = frozenset(
    {"Immune", "Endothelial", "Epithelial", "Neural", "Mesenchymal"}
)


@dataclass(frozen=True)
class FamilySpec:
    """A fixed multiple-testing family.

    ``familywise_alpha`` is divided equally over the declared members.  The
    allocation therefore depends only on the versioned registry, never on how
    many markers happen to be detected or expressed in a donor.
    """

    name: str
    familywise_alpha: float
    members: tuple[str, ...]

    @property
    def member_alpha(self) -> float:
        return self.familywise_alpha / len(self.members)


@dataclass(frozen=True)
class MarkerSpec:
    """One marker's interpretation and measurement policy."""

    name: str
    kind: str
    uses: tuple[str, ...]
    broad_compartment: str | None
    calibration_status: str
    reference: str | None
    families: tuple[str, ...]
    compartment: str
    spillover_policy: str
    mutually_exclusive: bool

    @property
    def active(self) -> bool:
        return self.calibration_status == "active"

    @property
    def subtype_only(self) -> bool:
        return "subtype" in self.uses and "broad" not in self.uses

    @property
    def measurement_column(self) -> str:
        """Corrected-compartment column selected by this specification."""

        return f"{self.compartment}__{self.name}"


@dataclass(frozen=True)
class MarkerRegistry:
    """Typed, validated view of ``marker_registry.json``."""

    schema_version: int
    registry_version: str
    method_version: str
    description: str
    allowed_compartments: tuple[str, ...]
    allowed_spillover_policies: tuple[str, ...]
    families: Mapping[str, FamilySpec]
    markers: tuple[MarkerSpec, ...]
    fingerprint: str

    @property
    def by_name(self) -> dict[str, MarkerSpec]:
        return {marker.name: marker for marker in self.markers}

    def marker(self, name: str) -> MarkerSpec:
        try:
            return self.by_name[name]
        except KeyError as exc:
            raise KeyError(f"marker is not registered: {name}") from exc

    @property
    def active_markers(self) -> tuple[MarkerSpec, ...]:
        return tuple(marker for marker in self.markers if marker.active)

    @property
    def broad_markers(self) -> tuple[MarkerSpec, ...]:
        return tuple(marker for marker in self.markers if "broad" in marker.uses)

    @property
    def type_markers(self) -> tuple[MarkerSpec, ...]:
        return tuple(marker for marker in self.markers if marker.kind == "type")

    @property
    def state_markers(self) -> tuple[MarkerSpec, ...]:
        return tuple(marker for marker in self.markers if marker.kind == "state")

    def markers_for_use(self, use: str) -> tuple[MarkerSpec, ...]:
        return tuple(marker for marker in self.markers if use in marker.uses)

    @property
    def reference_graph(self) -> dict[str, str]:
        """Directed active target -> mutually-exclusive reference graph.

        Reciprocal pairs (notably E-cadherin and CD31) are valid.  The graph is
        an audit representation, not an instruction to estimate one call from
        the other: calibration receives explicit reference-positive masks.
        """

        return {
            marker.name: marker.reference
            for marker in self.active_markers
            if marker.reference is not None
        }

    def targets_for_reference(self, reference: str) -> tuple[str, ...]:
        return tuple(
            marker.name
            for marker in self.active_markers
            if marker.reference == reference
        )

    def marker_alpha(self, marker: str | MarkerSpec) -> float:
        """Return the marker's fixed Bonferroni allocation.

        A marker can participate in more than one interpretation family (INS
        and GCG are both broad epithelial and endocrine subtype evidence).  In
        that case the strictest pre-declared allocation is used.
        """

        spec = self.marker(marker) if isinstance(marker, str) else marker
        if not spec.active:
            raise ValueError(f"{spec.name} has no active calibration")
        allocations = [self.families[name].member_alpha for name in spec.families]
        return min(allocations)


def _canonical_fingerprint(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _duplicates(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    duplicate: set[str] = set()
    for value in values:
        if value in seen:
            duplicate.add(value)
        seen.add(value)
    return sorted(duplicate)


def registry_from_dict(payload: Mapping[str, Any]) -> MarkerRegistry:
    """Parse and validate a registry mapping.

    The caller may use this for tests or for an explicitly versioned future
    registry.  Unknown top-level metadata is harmless, while every scientific
    field used by the calibrator is validated below.
    """

    try:
        family_rows = payload["families"]
        marker_rows = payload["markers"]
    except KeyError as exc:
        raise ValueError(f"registry is missing required key: {exc.args[0]}") from exc
    if not isinstance(family_rows, Mapping):
        raise ValueError("registry families must be an object")
    if not isinstance(marker_rows, list):
        raise ValueError("registry markers must be a list")

    families: dict[str, FamilySpec] = {}
    for name, row in family_rows.items():
        if not isinstance(row, Mapping):
            raise ValueError(f"family {name!r} must be an object")
        try:
            families[str(name)] = FamilySpec(
                name=str(name),
                familywise_alpha=float(row["familywise_alpha"]),
                members=tuple(str(value) for value in row["members"]),
            )
        except KeyError as exc:
            raise ValueError(
                f"family {name!r} is missing required key: {exc.args[0]}"
            ) from exc

    markers: list[MarkerSpec] = []
    required = {
        "name",
        "kind",
        "uses",
        "broad_compartment",
        "calibration_status",
        "reference",
        "families",
        "compartment",
        "spillover_policy",
        "mutually_exclusive",
    }
    for index, row in enumerate(marker_rows):
        if not isinstance(row, Mapping):
            raise ValueError(f"marker row {index} must be an object")
        missing = required.difference(row)
        if missing:
            raise ValueError(
                f"marker row {index} is missing required keys: {sorted(missing)}"
            )
        markers.append(
            MarkerSpec(
                name=str(row["name"]),
                kind=str(row["kind"]),
                uses=tuple(str(value) for value in row["uses"]),
                broad_compartment=(
                    None
                    if row["broad_compartment"] is None
                    else str(row["broad_compartment"])
                ),
                calibration_status=str(row["calibration_status"]),
                reference=(
                    None if row["reference"] is None else str(row["reference"])
                ),
                families=tuple(str(value) for value in row["families"]),
                compartment=str(row["compartment"]),
                spillover_policy=str(row["spillover_policy"]),
                mutually_exclusive=bool(row["mutually_exclusive"]),
            )
        )

    try:
        registry = MarkerRegistry(
            schema_version=int(payload["schema_version"]),
            registry_version=str(payload["registry_version"]),
            method_version=str(payload["method_version"]),
            description=str(payload.get("description", "")),
            allowed_compartments=tuple(
                str(value) for value in payload["allowed_compartments"]
            ),
            allowed_spillover_policies=tuple(
                str(value) for value in payload["allowed_spillover_policies"]
            ),
            families=families,
            markers=tuple(markers),
            fingerprint=_canonical_fingerprint(payload),
        )
    except KeyError as exc:
        raise ValueError(f"registry is missing required key: {exc.args[0]}") from exc
    validate_registry(registry)
    return registry


def validate_registry(registry: MarkerRegistry) -> None:
    """Raise ``ValueError`` for an internally contradictory registry."""

    errors: list[str] = []
    if registry.schema_version != 1:
        errors.append(f"unsupported schema_version {registry.schema_version!r}")
    if not registry.registry_version.strip():
        errors.append("registry_version must be non-empty")
    if not registry.method_version.strip():
        errors.append("method_version must be non-empty")
    if not registry.allowed_compartments:
        errors.append("allowed_compartments must be non-empty")
    if not registry.allowed_spillover_policies:
        errors.append("allowed_spillover_policies must be non-empty")

    names = [marker.name for marker in registry.markers]
    duplicated_names = _duplicates(names)
    if duplicated_names:
        errors.append(f"duplicate marker names: {duplicated_names}")
    name_set = set(names)

    for family in registry.families.values():
        if not (0.0 < family.familywise_alpha < 1.0):
            errors.append(
                f"family {family.name}: familywise_alpha must be between 0 and 1"
            )
        if not family.members:
            errors.append(f"family {family.name}: members must be non-empty")
        duplicate_members = _duplicates(family.members)
        if duplicate_members:
            errors.append(
                f"family {family.name}: duplicate members {duplicate_members}"
            )
        missing = sorted(set(family.members).difference(name_set))
        if missing:
            errors.append(f"family {family.name}: unknown members {missing}")

    for marker in registry.markers:
        label = f"marker {marker.name!r}"
        if not marker.name or marker.name.strip() != marker.name:
            errors.append(f"{label}: name must be non-empty and trimmed")
        if marker.kind not in ALLOWED_KINDS:
            errors.append(f"{label}: unknown kind {marker.kind!r}")
        if marker.calibration_status not in ALLOWED_CALIBRATION_STATUSES:
            errors.append(
                f"{label}: unknown calibration_status "
                f"{marker.calibration_status!r}"
            )
        if not marker.uses or _duplicates(marker.uses):
            errors.append(f"{label}: uses must be non-empty and unique")
        if marker.compartment not in registry.allowed_compartments:
            errors.append(
                f"{label}: unsupported compartment {marker.compartment!r}"
            )
        if marker.spillover_policy not in registry.allowed_spillover_policies:
            errors.append(
                f"{label}: unsupported spillover policy "
                f"{marker.spillover_policy!r}"
            )
        if marker.spillover_policy == "uncompensated" and marker.compartment != "Nucleus":
            errors.append(
                f"{label}: uncompensated per-cell markers must select Nucleus"
            )
        if marker.broad_compartment not in ALLOWED_BROAD_COMPARTMENTS | {None}:
            errors.append(
                f"{label}: unknown broad compartment {marker.broad_compartment!r}"
            )
        if "broad" in marker.uses and marker.broad_compartment is None:
            errors.append(f"{label}: broad use requires broad_compartment")
        if "broad" not in marker.uses and marker.broad_compartment is not None:
            # Deferred identity markers retain their intended broad compartment.
            if "deferred_type" not in marker.uses:
                errors.append(
                    f"{label}: broad_compartment without broad/deferred_type use"
                )

        if marker.active:
            if marker.reference is None:
                errors.append(f"{label}: active calibration requires a reference")
            elif marker.reference not in name_set:
                errors.append(
                    f"{label}: reference {marker.reference!r} is not registered"
                )
            elif marker.reference == marker.name:
                errors.append(f"{label}: marker cannot reference itself")
            if not marker.mutually_exclusive:
                errors.append(
                    f"{label}: active pair must be declared mutually exclusive"
                )
            if not marker.families:
                errors.append(
                    f"{label}: active calibration requires an error family"
                )
        else:
            if marker.reference is not None:
                errors.append(
                    f"{label}: inactive marker cannot carry an operative reference"
                )
            if marker.families:
                errors.append(
                    f"{label}: inactive marker cannot consume a testing allocation"
                )
            if marker.mutually_exclusive:
                errors.append(
                    f"{label}: inactive marker cannot claim a validated pair"
                )

        for family_name in marker.families:
            if family_name not in registry.families:
                errors.append(f"{label}: unknown family {family_name!r}")
            elif marker.name not in registry.families[family_name].members:
                errors.append(
                    f"{label}: family {family_name!r} does not list the marker"
                )

    by_name = registry.by_name
    for family in registry.families.values():
        for name in family.members:
            if name in by_name and family.name not in by_name[name].families:
                errors.append(
                    f"family {family.name}: member {name!r} does not point back "
                    "to the family"
                )
            if name in by_name and not by_name[name].active:
                errors.append(
                    f"family {family.name}: member {name!r} is not active"
                )

    # Locked decisions that prompted this registry.  Keeping them here makes a
    # stale edited JSON fail at load rather than silently reintroducing the
    # contradictory pair/gate paths this module replaces.
    cd11b = by_name.get("CD11b")
    if cd11b is None or cd11b.reference != "E_cadherin":
        errors.append("CD11b must use E_cadherin as its active reference")
    sst = by_name.get("SST")
    if (
        sst is None
        or sst.reference != "Pan_Cytokeratin"
        or not sst.subtype_only
        or not sst.active
    ):
        errors.append(
            "SST must be active subtype-only evidence with Pan_Cytokeratin "
            "as its reference"
        )

    if errors:
        raise ValueError("invalid marker registry:\n- " + "\n- ".join(errors))


def load_registry(path: str | Path | None = None) -> MarkerRegistry:
    """Load the default or explicitly supplied JSON registry."""

    source = DEFAULT_REGISTRY_PATH if path is None else Path(path)
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"marker registry not found: {source}") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON in marker registry {source}: {exc}") from exc
    if not isinstance(payload, Mapping):
        raise ValueError(f"marker registry {source} must contain a JSON object")
    return registry_from_dict(payload)


__all__ = [
    "ALLOWED_BROAD_COMPARTMENTS",
    "ALLOWED_CALIBRATION_STATUSES",
    "ALLOWED_KINDS",
    "DEFAULT_REGISTRY_PATH",
    "FamilySpec",
    "MarkerRegistry",
    "MarkerSpec",
    "load_registry",
    "registry_from_dict",
    "validate_registry",
]
