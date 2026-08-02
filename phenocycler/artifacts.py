"""Content-derived input and artifact contracts for the Phenocycler workflow.

The existing pipeline historically treated the presence of a donor directory as
proof that a stage was complete.  This module provides the stricter primitives a
workflow driver can use instead:

* immutable fingerprints for files, configuration and relevant source code;
* an explicit QuPath image contract (measurements, image, cell/nucleus geometry,
  pixel calibration, channel map and segmentation version);
* exact snapshots of donor-partitioned Parquet datasets, including schema and
  object-id-universe hashes; and
* stage/run manifests whose identifiers are derived from their content rather
  than their creation time.

The module intentionally has no dependency on ``pipeline.py``.  It can therefore
be introduced alongside the current pipeline and used by a replacement
orchestrator without making the old existence checks appear more trustworthy
than they are.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import subprocess
from dataclasses import asdict, dataclass, fields, is_dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence


MANIFEST_VERSION = 1
_HASH_CHUNK_BYTES = 8 * 1024 * 1024
_SAMPLE_BYTES = 1024 * 1024
_SAMPLE_COUNT = 3


class ContractError(ValueError):
    """Base class for a violated input or artifact contract."""


class GeometryContractError(ContractError):
    """A QuPath geometry export is absent, incomplete or malformed."""


class DuplicateObjectError(GeometryContractError):
    """A supposedly canonical object-id universe contains duplicates."""


class PartialArtifactError(ContractError):
    """The exact expected donor set is not complete."""


class SchemaContractError(ContractError):
    """Parquet files do not carry one exact, stable schema."""


class StaleArtifactError(ContractError):
    """A manifest no longer describes its inputs, code, config or outputs."""


def _normalise_json(value: Any) -> Any:
    """Convert common Python/scientific values to deterministic JSON values."""
    if is_dataclass(value):
        return _normalise_json(asdict(value))
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Mapping):
        return {
            str(key): _normalise_json(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_normalise_json(item) for item in value]
    if isinstance(value, (set, frozenset)):
        normalised = [_normalise_json(item) for item in value]
        return sorted(normalised, key=canonical_json)
    if hasattr(value, "item") and callable(value.item):
        try:
            return _normalise_json(value.item())
        except (TypeError, ValueError):
            pass
    if isinstance(value, float) and not math.isfinite(value):
        if math.isnan(value):
            return {"__nonfinite_float__": "nan"}
        return {"__nonfinite_float__": "inf" if value > 0 else "-inf"}
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def canonical_json(value: Any) -> str:
    """Return a deterministic, whitespace-free JSON representation."""
    return json.dumps(
        _normalise_json(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def sha256_json(value: Any) -> str:
    """Hash a value after deterministic JSON normalisation."""
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def _content_identity(value: Any) -> Any:
    """Remove validation-cache/timestamp fields from a content identifier."""
    normalised = _normalise_json(value)
    if isinstance(normalised, dict):
        return {
            key: _content_identity(item)
            for key, item in normalised.items()
            if key not in {"created_at_utc", "mtime_ns"}
        }
    if isinstance(normalised, list):
        return [_content_identity(item) for item in normalised]
    return normalised


def sha256_file(path: Path | str) -> str:
    """Hash all bytes in ``path`` (no mtime/size shortcut)."""
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(_HASH_CHUNK_BYTES):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_sampled_file(
    path: Path | str,
    *,
    sample_size_bytes: int = _SAMPLE_BYTES,
    sample_count: int = _SAMPLE_COUNT,
) -> str:
    """Hash deterministic windows across a large file.

    This is an explicitly weaker, opt-in fingerprint for environments where an
    initial full read of a multi-hundred-GB qptiff is infeasible. The file size,
    sample layout and offsets are included in the digest so the result cannot be
    confused with a full SHA-256.
    """
    path = Path(path)
    if sample_size_bytes <= 0 or sample_count <= 0:
        raise ContractError("sample size and count must be positive")
    size = path.stat().st_size
    width = min(sample_size_bytes, size)
    if size <= width or sample_count == 1:
        offsets = (0,)
    else:
        last = size - width
        offsets = tuple(
            sorted({round(index * last / (sample_count - 1)) for index in range(sample_count)})
        )
    digest = hashlib.sha256()
    digest.update(b"phenocycler-sampled-sha256-v1\0")
    digest.update(size.to_bytes(16, "big"))
    digest.update(sample_size_bytes.to_bytes(8, "big"))
    digest.update(sample_count.to_bytes(8, "big"))
    with path.open("rb") as handle:
        for offset in offsets:
            handle.seek(offset)
            chunk = handle.read(width)
            digest.update(offset.to_bytes(16, "big"))
            digest.update(len(chunk).to_bytes(8, "big"))
            digest.update(chunk)
    return digest.hexdigest()


def hash_identifiers(values: Iterable[Any]) -> str:
    """Hash an exact identifier set independent of row/file order.

    Length-prefixing prevents ambiguous concatenations (``["ab", "c"]`` versus
    ``["a", "bc"]``). Nulls and duplicates are contract violations.
    """
    ids: list[str] = []
    for value in values:
        if value is None:
            raise ContractError("object_id contains null")
        text = str(value)
        if not text:
            raise ContractError("object_id contains an empty value")
        ids.append(text)
    ids.sort()
    for previous, current in zip(ids, ids[1:]):
        if previous == current:
            raise DuplicateObjectError(f"duplicate object_id: {current}")
    digest = hashlib.sha256()
    for value in ids:
        encoded = value.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


def dataframe_content_sha256(frame, *, sort_by: Sequence[str] = ()) -> str:
    """Hash tabular values deterministically without relying on pickle internals."""
    missing = set(sort_by).difference(frame.columns)
    if missing:
        raise ContractError(f"cannot hash frame: missing sort columns {sorted(missing)}")
    ordered = frame.copy()
    if sort_by:
        ordered = ordered.sort_values(list(sort_by), kind="mergesort")
    rows = []
    for row in ordered.itertuples(index=False, name=None):
        values = []
        for value in row:
            try:
                missing_value = bool(value is None or value != value)
            except (TypeError, ValueError):
                missing_value = False
            values.append(None if missing_value else _normalise_json(value))
        rows.append(values)
    return sha256_json({"columns": list(map(str, ordered.columns)), "rows": rows})


@dataclass(frozen=True)
class FileFingerprint:
    """An exact path + byte-content fingerprint."""

    path: str
    size_bytes: int
    mtime_ns: int
    sha256: str
    hash_mode: str = "full"
    sample_size_bytes: int = 0
    sample_count: int = 0

    @classmethod
    def capture(
        cls,
        path: Path | str,
        *,
        hash_mode: str = "full",
        sample_size_bytes: int = _SAMPLE_BYTES,
        sample_count: int = _SAMPLE_COUNT,
    ) -> "FileFingerprint":
        resolved = Path(path).expanduser().resolve()
        if not resolved.is_file():
            raise ContractError(f"required file does not exist: {resolved}")
        if hash_mode not in {"full", "sampled"}:
            raise ContractError("hash_mode must be 'full' or 'sampled'")
        stat = resolved.stat()
        digest = (
            sha256_file(resolved)
            if hash_mode == "full"
            else sha256_sampled_file(
                resolved,
                sample_size_bytes=sample_size_bytes,
                sample_count=sample_count,
            )
        )
        return cls(
            path=str(resolved),
            size_bytes=stat.st_size,
            mtime_ns=stat.st_mtime_ns,
            sha256=digest,
            hash_mode=hash_mode,
            sample_size_bytes=sample_size_bytes if hash_mode == "sampled" else 0,
            sample_count=sample_count if hash_mode == "sampled" else 0,
        )

    def _current_digest(self) -> str:
        if self.hash_mode == "full":
            return sha256_file(self.path)
        if self.hash_mode == "sampled":
            return sha256_sampled_file(
                self.path,
                sample_size_bytes=self.sample_size_bytes,
                sample_count=self.sample_count,
            )
        raise ContractError(f"unsupported fingerprint mode {self.hash_mode!r}")

    def validate_current(self, *, mode: str = "fast") -> None:
        """Validate with ``fast`` (default), ``stat`` or ``content`` semantics.

        ``fast`` first compares cached path/size/mtime and returns immediately
        when they match. If only mtime changed, it rehashes once to distinguish a
        touched file from changed content. ``stat`` never reads file contents.
        ``content`` always recomputes the recorded full/sampled digest.
        """
        if mode not in {"fast", "stat", "content"}:
            raise ContractError("validation mode must be 'fast', 'stat' or 'content'")
        path = Path(self.path)
        if not path.is_file():
            raise StaleArtifactError(f"fingerprinted file is missing: {path}")
        stat = path.stat()
        size_matches = stat.st_size == self.size_bytes
        mtime_matches = stat.st_mtime_ns == self.mtime_ns
        if mode == "stat":
            if not size_matches or not mtime_matches:
                raise StaleArtifactError(f"file stat changed: {path}")
            return
        if mode == "fast" and size_matches and mtime_matches:
            return
        if not size_matches:
            raise StaleArtifactError(
                f"file size changed: {path} "
                f"(manifest={self.size_bytes}, current={stat.st_size})"
            )
        current_digest = self._current_digest()
        if current_digest != self.sha256:
            raise StaleArtifactError(
                f"file content changed: {path} "
                f"(manifest={self.sha256[:12]}, current={current_digest[:12]})"
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "FileFingerprint":
        return cls(
            path=str(value["path"]),
            size_bytes=int(value["size_bytes"]),
            mtime_ns=int(value["mtime_ns"]),
            sha256=str(value["sha256"]),
            hash_mode=str(value.get("hash_mode", "full")),
            sample_size_bytes=int(value.get("sample_size_bytes", 0)),
            sample_count=int(value.get("sample_count", 0)),
        )


def _tree_entries(path: Path) -> tuple[tuple[str, int, str], ...]:
    if path.is_file():
        return ((path.name, path.stat().st_size, sha256_file(path)),)
    if not path.is_dir():
        raise ContractError(f"input path does not exist: {path}")
    entries = []
    for item in sorted(candidate for candidate in path.rglob("*") if candidate.is_file()):
        relative = item.relative_to(path).as_posix()
        entries.append((relative, item.stat().st_size, sha256_file(item)))
    return tuple(entries)


def fingerprint_path(path: Path | str) -> str:
    """Return a content hash for a file or a recursive directory tree."""
    resolved = Path(path).expanduser().resolve()
    return sha256_json(_tree_entries(resolved))


@dataclass(frozen=True)
class ConfigSnapshot:
    """Canonical configuration text plus its content hash."""

    canonical: str
    sha256: str

    @classmethod
    def capture(cls, config: Any) -> "ConfigSnapshot":
        if is_dataclass(config):
            values = {field.name: getattr(config, field.name) for field in fields(config)}
        elif isinstance(config, Mapping):
            values = dict(config)
        elif hasattr(config, "__dict__"):
            values = {
                key: value
                for key, value in vars(config).items()
                if not key.startswith("_")
            }
        else:
            values = config
        canonical = canonical_json(values)
        return cls(
            canonical=canonical,
            sha256=hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ConfigSnapshot":
        snapshot = cls(canonical=str(value["canonical"]), sha256=str(value["sha256"]))
        if hashlib.sha256(snapshot.canonical.encode("utf-8")).hexdigest() != snapshot.sha256:
            raise ContractError("configuration snapshot hash is invalid")
        return snapshot


def _source_entries(
    repository: Path, source_paths: Sequence[Path | str]
) -> tuple[tuple[str, int, str], ...]:
    entries: list[tuple[str, int, str]] = []
    seen: set[Path] = set()
    for raw in source_paths:
        candidate = Path(raw)
        if not candidate.is_absolute():
            candidate = repository / candidate
        candidate = candidate.resolve()
        if candidate.is_file():
            files_to_hash = [candidate]
        elif candidate.is_dir():
            files_to_hash = sorted(
                item
                for item in candidate.rglob("*")
                if item.is_file()
                and ".git" not in item.parts
                and "__pycache__" not in item.parts
            )
        else:
            raise ContractError(f"source path does not exist: {candidate}")
        for item in files_to_hash:
            if item in seen:
                continue
            seen.add(item)
            try:
                relative = item.relative_to(repository).as_posix()
            except ValueError:
                relative = str(item)
            entries.append((relative, item.stat().st_size, sha256_file(item)))
    return tuple(sorted(entries))


def _git_output(repository: Path, *args: str) -> str | None:
    try:
        completed = subprocess.run(
            ["git", "-C", str(repository), *args],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None
    return completed.stdout


@dataclass(frozen=True)
class CodeMetadata:
    """Commit/worktree metadata plus a byte hash of relevant source files."""

    repository: str
    commit: str | None
    dirty: bool
    worktree_status_sha256: str | None
    source_paths: tuple[str, ...]
    source_sha256: str

    @classmethod
    def capture(
        cls,
        repository: Path | str,
        *,
        source_paths: Sequence[Path | str],
    ) -> "CodeMetadata":
        root = Path(repository).expanduser().resolve()
        if not root.is_dir():
            raise ContractError(f"code repository does not exist: {root}")
        if not source_paths:
            raise ContractError("source_paths must name the code that produced the artifact")
        commit_raw = _git_output(root, "rev-parse", "HEAD")
        status = _git_output(root, "status", "--porcelain=v1", "--untracked-files=all")
        entries = _source_entries(root, source_paths)
        return cls(
            repository=str(root),
            commit=commit_raw.strip() if commit_raw else None,
            dirty=bool(status),
            worktree_status_sha256=(
                hashlib.sha256(status.encode("utf-8")).hexdigest()
                if status is not None
                else None
            ),
            source_paths=tuple(str(path) for path in source_paths),
            source_sha256=sha256_json(entries),
        )

    def capture_current(self) -> "CodeMetadata":
        return CodeMetadata.capture(self.repository, source_paths=self.source_paths)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "CodeMetadata":
        return cls(
            repository=str(value["repository"]),
            commit=value.get("commit"),
            dirty=bool(value["dirty"]),
            worktree_status_sha256=value.get("worktree_status_sha256"),
            source_paths=tuple(map(str, value["source_paths"])),
            source_sha256=str(value["source_sha256"]),
        )


@dataclass(frozen=True)
class ChannelMapping:
    """One exact marker-to-qptiff-channel mapping."""

    marker: str
    channel_name: str
    channel_index: int

    def __post_init__(self) -> None:
        if not self.marker.strip():
            raise ContractError("channel-map marker cannot be empty")
        if not self.channel_name.strip():
            raise ContractError(f"{self.marker}: channel name cannot be empty")
        if self.channel_index < 0:
            raise ContractError(f"{self.marker}: channel index must be >= 0")

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ChannelMapping":
        return cls(
            marker=str(value["marker"]),
            channel_name=str(value["channel_name"]),
            channel_index=int(value["channel_index"]),
        )


@dataclass(frozen=True)
class GeoJSONSummary:
    feature_count: int
    primary_geometry_count: int
    nucleus_geometry_count: int
    object_id_sha256: str

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "GeoJSONSummary":
        return cls(
            feature_count=int(value["feature_count"]),
            primary_geometry_count=int(value["primary_geometry_count"]),
            nucleus_geometry_count=int(value["nucleus_geometry_count"]),
            object_id_sha256=str(value["object_id_sha256"]),
        )


def geometry_present(geometry: Any) -> bool:
    if not isinstance(geometry, Mapping) or not geometry.get("type"):
        return False
    if geometry.get("type") == "GeometryCollection":
        return bool(geometry.get("geometries"))
    return bool(geometry.get("coordinates"))


class _JSONStream:
    """Incrementally decode JSON values without retaining the complete file."""

    def __init__(self, path: Path, *, chunk_size: int = 1024 * 1024):
        self.path = path
        self.handle = path.open(encoding="utf-8")
        self.chunk_size = chunk_size
        self.buffer = ""
        self.position = 0
        self.eof = False
        self.decoder = json.JSONDecoder()

    def close(self) -> None:
        self.handle.close()

    def _read_more(self) -> bool:
        if self.eof:
            return False
        if self.position:
            self.buffer = self.buffer[self.position :]
            self.position = 0
        chunk = self.handle.read(self.chunk_size)
        if not chunk:
            self.eof = True
            return False
        self.buffer += chunk
        return True

    def _skip_whitespace(self) -> None:
        while True:
            while (
                self.position < len(self.buffer)
                and self.buffer[self.position].isspace()
            ):
                self.position += 1
            if self.position < len(self.buffer) or not self._read_more():
                return

    def peek(self) -> str:
        self._skip_whitespace()
        if self.position >= len(self.buffer):
            raise json.JSONDecodeError(
                "unexpected end of JSON", self.buffer, self.position
            )
        return self.buffer[self.position]

    def expect(self, token: str) -> None:
        if self.peek() != token:
            raise json.JSONDecodeError(
                f"expected {token!r}", self.buffer, self.position
            )
        self.position += 1

    def decode(self) -> Any:
        self._skip_whitespace()
        while True:
            try:
                value, end = self.decoder.raw_decode(self.buffer, self.position)
            except json.JSONDecodeError:
                if not self._read_more():
                    raise
            else:
                self.position = end
                return value

    def require_end(self) -> None:
        self._skip_whitespace()
        if self.position < len(self.buffer) or self._read_more():
            self._skip_whitespace()
            if self.position < len(self.buffer):
                raise json.JSONDecodeError(
                    "extra data", self.buffer, self.position
                )


def iter_qupath_geojson_features(path: Path | str) -> Iterator[Mapping[str, Any]]:
    """Stream one QuPath FeatureCollection feature at a time.

    Cell exports in this cohort are several gigabytes each.  Using
    :func:`json.load` expands the full coordinate payload in memory even though
    manifest inspection and rasterization need only one feature at a time.
    """

    path = Path(path)
    if not path.is_file():
        raise GeometryContractError(f"GeoJSON is missing: {path}")
    stream = _JSONStream(path)
    top_level_type: Any = None
    found_features = False
    try:
        stream.expect("{")
        first_member = True
        while stream.peek() != "}":
            if not first_member:
                stream.expect(",")
            key = stream.decode()
            if not isinstance(key, str):
                raise GeometryContractError(
                    f"{path}: top-level GeoJSON key is not a string"
                )
            stream.expect(":")
            if key == "type":
                top_level_type = stream.decode()
            elif key == "features":
                if found_features:
                    raise GeometryContractError(
                        f"{path}: duplicate top-level features member"
                    )
                found_features = True
                stream.expect("[")
                first_feature = True
                while stream.peek() != "]":
                    if not first_feature:
                        stream.expect(",")
                    feature = stream.decode()
                    if not isinstance(feature, Mapping):
                        raise GeometryContractError(
                            f"{path}: feature is not an object"
                        )
                    yield feature
                    first_feature = False
                stream.expect("]")
            else:
                stream.decode()
            first_member = False
        stream.expect("}")
        stream.require_end()
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise GeometryContractError(f"cannot read GeoJSON {path}: {exc}") from exc
    finally:
        stream.close()
    if top_level_type != "FeatureCollection" or not found_features:
        raise GeometryContractError(f"{path}: expected a GeoJSON FeatureCollection")


def inspect_qupath_geojson(
    path: Path | str,
    *,
    require_primary_geometry: bool = True,
    require_nucleus_geometry: bool = False,
) -> tuple[GeoJSONSummary, tuple[str, ...]]:
    """Validate a QuPath FeatureCollection and return its exact object universe."""
    path = Path(path)
    object_ids: list[str] = []
    primary_count = 0
    nucleus_count = 0
    for index, feature in enumerate(iter_qupath_geojson_features(path)):
        object_id = feature.get("id")
        if object_id is None or not str(object_id):
            raise GeometryContractError(f"{path}: feature {index} has no object id")
        object_ids.append(str(object_id))
        primary = geometry_present(feature.get("geometry"))
        nucleus = geometry_present(feature.get("nucleusGeometry"))
        primary_count += int(primary)
        nucleus_count += int(nucleus)
        if require_primary_geometry and not primary:
            raise GeometryContractError(
                f"{path}: object {object_id} has no usable cell geometry"
            )
        if require_nucleus_geometry and not nucleus:
            raise GeometryContractError(
                f"{path}: object {object_id} has no usable nucleusGeometry"
            )
    object_hash = hash_identifiers(object_ids)
    return (
        GeoJSONSummary(
            feature_count=len(object_ids),
            primary_geometry_count=primary_count,
            nucleus_geometry_count=nucleus_count,
            object_id_sha256=object_hash,
        ),
        tuple(object_ids),
    )


@dataclass(frozen=True)
class QuPathImageManifest:
    """The complete immutable contract for one QuPath image export."""

    manifest_version: int
    donor_id: str
    image_id: str
    measurement_csv: FileFingerprint
    qptiff: FileFingerprint
    cell_geojson: FileFingerprint
    nucleus_geojson: FileFingerprint
    pixel_size_um_x: float
    pixel_size_um_y: float
    panel_id: str
    channel_map: tuple[ChannelMapping, ...]
    segmentation_version: str
    measurement_columns: tuple[str, ...]
    cell_geometry: GeoJSONSummary
    nucleus_geometry: GeoJSONSummary
    combined_geometry_file: bool
    content_id: str

    @classmethod
    def create(
        cls,
        *,
        donor_id: str,
        image_id: str,
        measurement_csv: Path | str,
        qptiff: Path | str,
        cell_geojson: Path | str,
        nucleus_geojson: Path | str | None,
        pixel_size_um_x: float,
        pixel_size_um_y: float,
        panel_id: str,
        channel_map: Sequence[ChannelMapping],
        segmentation_version: str,
        qptiff_hash_mode: str = "full",
        qptiff_sample_size_bytes: int = _SAMPLE_BYTES,
        qptiff_sample_count: int = _SAMPLE_COUNT,
        measurement_fingerprint: FileFingerprint | None = None,
        qptiff_fingerprint: FileFingerprint | None = None,
    ) -> "QuPathImageManifest":
        donor_id = str(donor_id).strip()
        image_id = str(image_id).strip()
        if not donor_id or not image_id:
            raise ContractError("donor_id and image_id are required")
        if pixel_size_um_x <= 0 or pixel_size_um_y <= 0:
            raise ContractError("pixel calibration must be finite and positive")
        if not math.isfinite(pixel_size_um_x) or not math.isfinite(pixel_size_um_y):
            raise ContractError("pixel calibration must be finite and positive")
        if not panel_id.strip() or not segmentation_version.strip():
            raise ContractError("panel_id and segmentation_version are required")
        channels = tuple(channel_map)
        if not channels:
            raise ContractError("channel_map cannot be empty")
        markers = [channel.marker for channel in channels]
        names = [channel.channel_name for channel in channels]
        indices = [channel.channel_index for channel in channels]
        if len(set(markers)) != len(markers):
            raise ContractError("channel_map contains duplicate markers")
        if len(set(names)) != len(names) or len(set(indices)) != len(indices):
            raise ContractError("channel_map maps more than one marker to a channel")

        measurement_path = Path(measurement_csv).expanduser().resolve()
        if measurement_fingerprint is None:
            measurement_fp = FileFingerprint.capture(measurement_path)
        else:
            measurement_fp = measurement_fingerprint
            if Path(measurement_fp.path) != measurement_path:
                raise ContractError(
                    "precomputed measurement fingerprint path differs from "
                    f"measurement_csv: {measurement_fp.path} != {measurement_path}"
                )
            measurement_fp.validate_current(mode="fast")
        with Path(measurement_fp.path).open(newline="", encoding="utf-8-sig") as handle:
            try:
                measurement_columns = tuple(next(csv.reader(handle)))
            except StopIteration as exc:
                raise ContractError("QuPath measurement CSV is empty") from exc
        for required in ("Image", "Object ID"):
            if required not in measurement_columns:
                raise ContractError(
                    f"QuPath measurement CSV is missing required column {required!r}"
                )

        cell_path = Path(cell_geojson).expanduser().resolve()
        nucleus_path = (
            Path(nucleus_geojson).expanduser().resolve()
            if nucleus_geojson is not None
            else cell_path
        )
        combined = cell_path == nucleus_path
        if combined:
            cell_summary, cell_ids = inspect_qupath_geojson(
                cell_path,
                require_primary_geometry=True,
            )
            if cell_summary.nucleus_geometry_count == 0:
                raise GeometryContractError(
                    f"{cell_path}: no feature has a usable nucleusGeometry"
                )
            nucleus_summary = GeoJSONSummary(
                feature_count=cell_summary.feature_count,
                primary_geometry_count=cell_summary.nucleus_geometry_count,
                nucleus_geometry_count=cell_summary.nucleus_geometry_count,
                object_id_sha256=cell_summary.object_id_sha256,
            )
        else:
            cell_summary, cell_ids = inspect_qupath_geojson(
                cell_path, require_primary_geometry=True
            )
            nucleus_summary, nucleus_ids = inspect_qupath_geojson(
                nucleus_path, require_primary_geometry=True
            )
            if set(cell_ids) != set(nucleus_ids):
                raise GeometryContractError(
                    "cell and nucleus GeoJSON object-id universes differ"
                )
        draft = cls(
            manifest_version=MANIFEST_VERSION,
            donor_id=donor_id,
            image_id=image_id,
            measurement_csv=measurement_fp,
            qptiff=(
                FileFingerprint.capture(
                    qptiff,
                    hash_mode=qptiff_hash_mode,
                    sample_size_bytes=qptiff_sample_size_bytes,
                    sample_count=qptiff_sample_count,
                )
                if qptiff_fingerprint is None
                else qptiff_fingerprint
            ),
            cell_geojson=FileFingerprint.capture(cell_path),
            nucleus_geojson=FileFingerprint.capture(nucleus_path),
            pixel_size_um_x=float(pixel_size_um_x),
            pixel_size_um_y=float(pixel_size_um_y),
            panel_id=panel_id.strip(),
            channel_map=tuple(sorted(channels, key=lambda item: item.channel_index)),
            segmentation_version=segmentation_version.strip(),
            measurement_columns=measurement_columns,
            cell_geometry=cell_summary,
            nucleus_geometry=nucleus_summary,
            combined_geometry_file=combined,
            content_id="",
        )
        if Path(draft.qptiff.path) != Path(qptiff).expanduser().resolve():
            raise ContractError(
                "precomputed qptiff fingerprint path differs from qptiff: "
                f"{draft.qptiff.path} != {Path(qptiff).expanduser().resolve()}"
            )
        draft.qptiff.validate_current(mode="fast")
        return replace(draft, content_id=sha256_json(draft._identity_payload()))

    def _identity_payload(self) -> dict[str, Any]:
        value = self.to_dict()
        value.pop("content_id", None)
        return _content_identity(value)

    def to_dict(self) -> dict[str, Any]:
        return _normalise_json(asdict(self))

    def verify_identity(self) -> None:
        if self.manifest_version != MANIFEST_VERSION:
            raise ContractError(
                f"unsupported QuPath manifest version {self.manifest_version}"
            )
        expected = sha256_json(self._identity_payload())
        if expected != self.content_id:
            raise ContractError("QuPath manifest content_id is invalid")

    def validate_current(self, *, mode: str = "fast") -> None:
        """Validate source files.

        ``mode="fast"`` is intended for routine status: unchanged qptiffs are
        accepted from cached path/size/mtime without reading their content.
        ``mode="content"`` performs a full (or explicitly sampled) rehash and
        repeats the geometry inspection.
        """
        self.verify_identity()
        for fingerprint in (
            self.measurement_csv,
            self.qptiff,
            self.cell_geojson,
            self.nucleus_geojson,
        ):
            fingerprint.validate_current(mode=mode)
        if mode != "content":
            # The unchanged content digest already covers the stored geometry
            # summary; reparsing a multi-GB FeatureCollection would add no
            # information to the fast status path.
            return
        if self.combined_geometry_file:
            cell_summary, _ = inspect_qupath_geojson(
                self.cell_geojson.path,
                require_primary_geometry=True,
            )
            if cell_summary.nucleus_geometry_count == 0:
                raise GeometryContractError(
                    f"{self.cell_geojson.path}: no feature has a usable "
                    "nucleusGeometry"
                )
            if cell_summary != self.cell_geometry:
                raise StaleArtifactError("combined cell/nucleus geometry summary changed")
        else:
            cell_summary, cell_ids = inspect_qupath_geojson(
                self.cell_geojson.path, require_primary_geometry=True
            )
            nucleus_summary, nucleus_ids = inspect_qupath_geojson(
                self.nucleus_geojson.path, require_primary_geometry=True
            )
            if set(cell_ids) != set(nucleus_ids):
                raise GeometryContractError(
                    "cell and nucleus GeoJSON object-id universes differ"
                )
            if cell_summary != self.cell_geometry or nucleus_summary != self.nucleus_geometry:
                raise StaleArtifactError("cell or nucleus geometry summary changed")

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "QuPathImageManifest":
        result = cls(
            manifest_version=int(value["manifest_version"]),
            donor_id=str(value["donor_id"]),
            image_id=str(value["image_id"]),
            measurement_csv=FileFingerprint.from_dict(value["measurement_csv"]),
            qptiff=FileFingerprint.from_dict(value["qptiff"]),
            cell_geojson=FileFingerprint.from_dict(value["cell_geojson"]),
            nucleus_geojson=FileFingerprint.from_dict(value["nucleus_geojson"]),
            pixel_size_um_x=float(value["pixel_size_um_x"]),
            pixel_size_um_y=float(value["pixel_size_um_y"]),
            panel_id=str(value["panel_id"]),
            channel_map=tuple(
                ChannelMapping.from_dict(item) for item in value["channel_map"]
            ),
            segmentation_version=str(value["segmentation_version"]),
            measurement_columns=tuple(map(str, value["measurement_columns"])),
            cell_geometry=GeoJSONSummary.from_dict(value["cell_geometry"]),
            nucleus_geometry=GeoJSONSummary.from_dict(value["nucleus_geometry"]),
            combined_geometry_file=bool(value["combined_geometry_file"]),
            content_id=str(value["content_id"]),
        )
        result.verify_identity()
        return result

    @classmethod
    def read_json(cls, path: Path | str) -> "QuPathImageManifest":
        with Path(path).open(encoding="utf-8") as handle:
            return cls.from_dict(json.load(handle))


@dataclass(frozen=True)
class PanelChannelMap:
    """The one exact channel mapping permitted for a named acquisition panel."""

    panel_id: str
    channels: tuple[ChannelMapping, ...]
    content_id: str

    @classmethod
    def create(
        cls, panel_id: str, channels: Sequence[ChannelMapping]
    ) -> "PanelChannelMap":
        panel_id = str(panel_id).strip()
        if not panel_id:
            raise ContractError("panel_id cannot be empty")
        ordered = tuple(sorted(channels, key=lambda item: item.channel_index))
        if not ordered:
            raise ContractError(f"panel {panel_id!r} has no channels")
        markers = [channel.marker for channel in ordered]
        names = [channel.channel_name for channel in ordered]
        indices = [channel.channel_index for channel in ordered]
        if len(set(markers)) != len(markers):
            raise ContractError(f"panel {panel_id!r} contains duplicate markers")
        if len(set(names)) != len(names) or len(set(indices)) != len(indices):
            raise ContractError(
                f"panel {panel_id!r} maps multiple markers to one channel"
            )
        payload = {
            "panel_id": panel_id,
            "channels": [asdict(channel) for channel in ordered],
        }
        return cls(
            panel_id=panel_id,
            channels=ordered,
            content_id=sha256_json(payload),
        )

    def verify_identity(self) -> None:
        expected = PanelChannelMap.create(self.panel_id, self.channels)
        if expected.content_id != self.content_id:
            raise ContractError(f"panel {self.panel_id!r} content_id is invalid")

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "PanelChannelMap":
        result = cls(
            panel_id=str(value["panel_id"]),
            channels=tuple(
                ChannelMapping.from_dict(item) for item in value["channels"]
            ),
            content_id=str(value["content_id"]),
        )
        result.verify_identity()
        return result


def _cohort_panel_maps(
    images: Sequence[QuPathImageManifest],
) -> tuple[PanelChannelMap, ...]:
    by_panel: dict[str, tuple[ChannelMapping, ...]] = {}
    for image in images:
        mapping = tuple(
            sorted(image.channel_map, key=lambda item: item.channel_index)
        )
        previous = by_panel.get(image.panel_id)
        if previous is not None and previous != mapping:
            raise ContractError(
                f"panel {image.panel_id!r} has contradictory channel mappings "
                f"between QuPath images"
            )
        by_panel[image.panel_id] = mapping
    return tuple(
        PanelChannelMap.create(panel_id, by_panel[panel_id])
        for panel_id in sorted(by_panel)
    )


def _validate_cohort_images(
    images: Sequence[QuPathImageManifest],
    expected_donors: Sequence[str],
) -> tuple[tuple[QuPathImageManifest, ...], tuple[str, ...]]:
    records = tuple(
        sorted(images, key=lambda item: (item.donor_id, item.image_id))
    )
    if not records:
        raise ContractError("a QuPath cohort manifest requires at least one image")
    for image in records:
        image.verify_identity()

    donor_ids = [image.donor_id for image in records]
    image_ids = [image.image_id for image in records]
    duplicate_donors = sorted(
        donor for donor in set(donor_ids) if donor_ids.count(donor) > 1
    )
    duplicate_images = sorted(
        image for image in set(image_ids) if image_ids.count(image) > 1
    )
    if duplicate_donors:
        raise ContractError(
            "QuPath cohort requires exactly one image per donor; duplicate donors="
            f"{duplicate_donors}"
        )
    if duplicate_images:
        raise ContractError(
            f"QuPath cohort contains duplicate image_id values={duplicate_images}"
        )

    expected = tuple(sorted(map(str, expected_donors)))
    if not expected or len(set(expected)) != len(expected):
        raise ContractError("expected_donors must be a non-empty unique set")
    actual = set(donor_ids)
    expected_set = set(expected)
    missing = sorted(expected_set - actual)
    extra = sorted(actual - expected_set)
    if missing or extra:
        raise PartialArtifactError(
            f"QuPath cohort donor set differs: missing={missing}, extra={extra}"
        )
    return records, expected


@dataclass(frozen=True)
class QuPathCohortManifest:
    """Exact cohort-level collection of one QuPath image per donor.

    Image manifests are embedded rather than linked indirectly, so the cohort
    content identifier owns every source digest and geometry contract. Images
    sharing a ``panel_id`` must carry byte-for-byte equivalent marker, channel
    name and channel-index mappings. Different named panels may intentionally
    differ.
    """

    manifest_version: int
    cohort_name: str
    expected_donors: tuple[str, ...]
    images: tuple[QuPathImageManifest, ...]
    panels: tuple[PanelChannelMap, ...]
    created_at_utc: str
    content_id: str

    @classmethod
    def create(
        cls,
        *,
        cohort_name: str,
        expected_donors: Sequence[str],
        images: Sequence[QuPathImageManifest],
        validate_sources: bool = True,
        validation_mode: str = "fast",
    ) -> "QuPathCohortManifest":
        cohort_name = str(cohort_name).strip()
        if not cohort_name:
            raise ContractError("cohort_name is required")
        records, expected = _validate_cohort_images(images, expected_donors)
        panels = _cohort_panel_maps(records)
        if validate_sources:
            for image in records:
                image.validate_current(mode=validation_mode)
        draft = cls(
            manifest_version=MANIFEST_VERSION,
            cohort_name=cohort_name,
            expected_donors=expected,
            images=records,
            panels=panels,
            created_at_utc=datetime.now(timezone.utc).isoformat(),
            content_id="",
        )
        return replace(draft, content_id=sha256_json(draft._identity_payload()))

    @property
    def donor_ids(self) -> tuple[str, ...]:
        return tuple(image.donor_id for image in self.images)

    @property
    def image_ids(self) -> tuple[str, ...]:
        return tuple(image.image_id for image in self.images)

    def image_to_donor(self) -> dict[str, str]:
        """Return a fresh ingest mapping without reading any source file."""
        return image_id_to_donor_id(self)

    def _identity_payload(self) -> dict[str, Any]:
        value = self.to_dict()
        value.pop("content_id", None)
        return _content_identity(value)

    def to_dict(self) -> dict[str, Any]:
        return _normalise_json(asdict(self))

    def verify_identity(self) -> None:
        if self.manifest_version != MANIFEST_VERSION:
            raise ContractError(
                f"unsupported QuPath cohort manifest version {self.manifest_version}"
            )
        records, expected = _validate_cohort_images(
            self.images, self.expected_donors
        )
        if records != self.images or expected != self.expected_donors:
            raise ContractError("QuPath cohort records are not canonically ordered")
        derived_panels = _cohort_panel_maps(self.images)
        for panel in self.panels:
            panel.verify_identity()
        if self.panels != derived_panels:
            raise ContractError(
                "QuPath cohort panel registry does not match its image manifests"
            )
        if sha256_json(self._identity_payload()) != self.content_id:
            raise ContractError("QuPath cohort manifest content_id is invalid")

    def validate_current(self, *, mode: str = "fast") -> None:
        """Validate every embedded image source using fast or content mode."""
        self.verify_identity()
        for image in self.images:
            image.validate_current(mode=mode)

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "QuPathCohortManifest":
        result = cls(
            manifest_version=int(value["manifest_version"]),
            cohort_name=str(value["cohort_name"]),
            expected_donors=tuple(map(str, value["expected_donors"])),
            images=tuple(
                QuPathImageManifest.from_dict(item) for item in value["images"]
            ),
            panels=tuple(
                PanelChannelMap.from_dict(item) for item in value["panels"]
            ),
            created_at_utc=str(value["created_at_utc"]),
            content_id=str(value["content_id"]),
        )
        result.verify_identity()
        return result

    @classmethod
    def read_json(cls, path: Path | str) -> "QuPathCohortManifest":
        with Path(path).open(encoding="utf-8") as handle:
            return cls.from_dict(json.load(handle))


def image_id_to_donor_id(
    cohort_or_images: QuPathCohortManifest | Sequence[QuPathImageManifest],
) -> dict[str, str]:
    """Purely derive the exact QuPath ``Image`` value -> donor mapping.

    No path is opened and no global state is consulted. This is the intended
    ingest replacement for extracting the first digit run from an image name.
    """
    images = (
        cohort_or_images.images
        if isinstance(cohort_or_images, QuPathCohortManifest)
        else tuple(cohort_or_images)
    )
    result: dict[str, str] = {}
    for image in images:
        if image.image_id in result:
            raise ContractError(
                f"duplicate image_id cannot map to one donor: {image.image_id!r}"
            )
        result[image.image_id] = image.donor_id
    return {image_id: result[image_id] for image_id in sorted(result)}


@dataclass(frozen=True)
class SchemaField:
    name: str
    type: str
    nullable: bool

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "SchemaField":
        return cls(
            name=str(value["name"]),
            type=str(value["type"]),
            nullable=bool(value["nullable"]),
        )


@dataclass(frozen=True)
class DonorDatasetSnapshot:
    donor_id: str
    files: tuple[FileFingerprint, ...]
    row_count: int
    schema: tuple[SchemaField, ...]
    schema_sha256: str
    object_id_sha256: str
    content_sha256: str

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "DonorDatasetSnapshot":
        return cls(
            donor_id=str(value["donor_id"]),
            files=tuple(FileFingerprint.from_dict(item) for item in value["files"]),
            row_count=int(value["row_count"]),
            schema=tuple(SchemaField.from_dict(item) for item in value["schema"]),
            schema_sha256=str(value["schema_sha256"]),
            object_id_sha256=str(value["object_id_sha256"]),
            content_sha256=str(value["content_sha256"]),
        )


@dataclass(frozen=True)
class DatasetSnapshot:
    root: str
    donors: tuple[DonorDatasetSnapshot, ...]
    total_rows: int
    schema: tuple[SchemaField, ...]
    schema_sha256: str
    object_id_sha256: str
    content_sha256: str

    @property
    def donor_ids(self) -> tuple[str, ...]:
        return tuple(donor.donor_id for donor in self.donors)

    def to_dict(self) -> dict[str, Any]:
        return _normalise_json(asdict(self))

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "DatasetSnapshot":
        return cls(
            root=str(value["root"]),
            donors=tuple(
                DonorDatasetSnapshot.from_dict(item) for item in value["donors"]
            ),
            total_rows=int(value["total_rows"]),
            schema=tuple(SchemaField.from_dict(item) for item in value["schema"]),
            schema_sha256=str(value["schema_sha256"]),
            object_id_sha256=str(value["object_id_sha256"]),
            content_sha256=str(value["content_sha256"]),
        )


def _schema_fields(schema) -> tuple[SchemaField, ...]:
    return tuple(
        SchemaField(name=field.name, type=str(field.type), nullable=field.nullable)
        for field in schema
    )


def snapshot_parquet_dataset(
    root: Path | str,
    *,
    expected_donors: Sequence[str] | None = None,
    object_id_column: str = "object_id",
    partition_prefix: str = "donor_id=",
) -> DatasetSnapshot:
    """Create an exact snapshot of a donor-partitioned Parquet dataset.

    If ``expected_donors`` is supplied, missing *and extra* donor partitions are
    fatal. Every Parquet file must have the same Arrow schema, each donor's
    object IDs must be unique, and an in-file ``donor_id`` column (when present)
    must agree with its partition.
    """
    import pyarrow.parquet as pq

    root = Path(root).expanduser().resolve()
    if not root.is_dir():
        raise PartialArtifactError(f"dataset root is missing: {root}")
    partitions = {
        part.name[len(partition_prefix) :]: part
        for part in root.iterdir()
        if part.is_dir() and part.name.startswith(partition_prefix)
    }
    actual = tuple(sorted(partitions))
    if expected_donors is not None:
        expected = tuple(sorted(map(str, expected_donors)))
        missing = sorted(set(expected).difference(actual))
        extra = sorted(set(actual).difference(expected))
        if missing or extra:
            raise PartialArtifactError(
                f"dataset donor set differs: missing={missing}, extra={extra}"
            )
    if not actual:
        raise PartialArtifactError(f"dataset has no donor partitions: {root}")

    donor_snapshots: list[DonorDatasetSnapshot] = []
    cohort_schema = None
    cohort_schema_hash = None
    cohort_object_parts = []
    for donor_id in actual:
        partition = partitions[donor_id]
        files = sorted(partition.glob("*.parquet"))
        if not files:
            raise PartialArtifactError(f"donor {donor_id} has no Parquet files")
        ids: list[Any] = []
        file_fingerprints = []
        row_count = 0
        donor_schema = None
        donor_schema_hash = None
        file_content_entries = []
        for path in files:
            parquet = pq.ParquetFile(path)
            schema = parquet.schema_arrow
            schema_hash = hashlib.sha256(schema.serialize().to_pybytes()).hexdigest()
            if object_id_column not in schema.names:
                raise SchemaContractError(
                    f"{path}: schema has no {object_id_column!r} column"
                )
            if donor_schema is None:
                donor_schema = schema
                donor_schema_hash = schema_hash
            elif not donor_schema.equals(schema, check_metadata=True):
                raise SchemaContractError(
                    f"donor {donor_id} contains inconsistent Parquet schemas"
                )
            if cohort_schema is None:
                cohort_schema = schema
                cohort_schema_hash = schema_hash
            elif not cohort_schema.equals(schema, check_metadata=True):
                raise SchemaContractError(
                    f"donor {donor_id} schema differs from the cohort schema"
                )
            columns = [object_id_column]
            if "donor_id" in schema.names:
                columns.append("donor_id")
            table = parquet.read(columns=columns)
            current_ids = table.column(object_id_column).to_pylist()
            ids.extend(current_ids)
            row_count += table.num_rows
            if "donor_id" in columns:
                stored = {
                    str(value)
                    for value in table.column("donor_id").to_pylist()
                    if value is not None
                }
                if stored != {donor_id}:
                    raise ContractError(
                        f"{path}: stored donor_id values {sorted(stored)} "
                        f"do not equal partition {donor_id}"
                    )
            fingerprint = FileFingerprint.capture(path)
            file_fingerprints.append(fingerprint)
            file_content_entries.append(
                (
                    path.relative_to(partition).as_posix(),
                    fingerprint.size_bytes,
                    fingerprint.sha256,
                )
            )
        object_hash = hash_identifiers(ids)
        donor_payload = {
            "donor_id": donor_id,
            "files": file_content_entries,
            "rows": row_count,
            "schema": donor_schema_hash,
            "objects": object_hash,
        }
        donor_content_hash = sha256_json(donor_payload)
        donor_snapshots.append(
            DonorDatasetSnapshot(
                donor_id=donor_id,
                files=tuple(file_fingerprints),
                row_count=row_count,
                schema=_schema_fields(donor_schema),
                schema_sha256=str(donor_schema_hash),
                object_id_sha256=object_hash,
                content_sha256=donor_content_hash,
            )
        )
        cohort_object_parts.append((donor_id, row_count, object_hash))

    total_rows = sum(donor.row_count for donor in donor_snapshots)
    cohort_object_hash = sha256_json(cohort_object_parts)
    content_hash = sha256_json(
        {
            "donors": [
                (donor.donor_id, donor.content_sha256) for donor in donor_snapshots
            ],
            "rows": total_rows,
            "schema": cohort_schema_hash,
            "objects": cohort_object_hash,
        }
    )
    return DatasetSnapshot(
        root=str(root),
        donors=tuple(donor_snapshots),
        total_rows=total_rows,
        schema=_schema_fields(cohort_schema),
        schema_sha256=str(cohort_schema_hash),
        object_id_sha256=cohort_object_hash,
        content_sha256=content_hash,
    )


def validate_dataset_snapshot_current(
    snapshot: DatasetSnapshot, *, mode: str = "fast"
) -> None:
    """Validate a dataset snapshot without re-reading unchanged Parquet content.

    Fast validation checks the exact donor/file inventory and every cached
    path/size/mtime. A file with unchanged metadata is not opened. If metadata
    changed, ``FileFingerprint`` confirms its content digest before accepting
    it. ``mode="content"`` additionally rebuilds schema/object/content hashes.
    """
    if mode not in {"fast", "stat", "content"}:
        raise ContractError("validation mode must be 'fast', 'stat' or 'content'")
    root = Path(snapshot.root)
    if not root.is_dir():
        raise PartialArtifactError(f"dataset root is missing: {root}")
    actual_donors = {
        path.name.split("=", 1)[1]
        for path in root.glob("donor_id=*")
        if path.is_dir()
    }
    expected_donors = set(snapshot.donor_ids)
    if actual_donors != expected_donors:
        raise PartialArtifactError(
            "dataset donor set differs: "
            f"missing={sorted(expected_donors - actual_donors)}, "
            f"extra={sorted(actual_donors - expected_donors)}"
        )
    for donor in snapshot.donors:
        partition = root / f"donor_id={donor.donor_id}"
        actual_files = {
            str(path.resolve()) for path in partition.glob("*.parquet")
        }
        recorded_files = {fingerprint.path for fingerprint in donor.files}
        if actual_files != recorded_files:
            raise PartialArtifactError(
                f"donor {donor.donor_id} Parquet file inventory changed"
            )
        for fingerprint in donor.files:
            fingerprint.validate_current(mode=mode)
    if mode == "content":
        current = snapshot_parquet_dataset(
            root, expected_donors=snapshot.donor_ids
        )
        if current.content_sha256 != snapshot.content_sha256:
            raise StaleArtifactError("dataset content/schema/object universe changed")


@dataclass(frozen=True)
class InputArtifact:
    """A content-addressed stage input."""

    name: str
    kind: str
    path: str
    content_sha256: str
    files: tuple[FileFingerprint, ...]

    @classmethod
    def from_path(
        cls,
        name: str,
        path: Path | str,
        *,
        hash_mode: str = "full",
    ) -> "InputArtifact":
        resolved = Path(path).expanduser().resolve()
        if not resolved.exists():
            raise ContractError(f"input artifact is missing: {resolved}")
        paths = (
            [resolved]
            if resolved.is_file()
            else sorted(item for item in resolved.rglob("*") if item.is_file())
        )
        fingerprints = tuple(
            FileFingerprint.capture(item, hash_mode=hash_mode) for item in paths
        )
        entries = [
            (
                (
                    Path(fingerprint.path).name
                    if resolved.is_file()
                    else Path(fingerprint.path).relative_to(resolved).as_posix()
                ),
                fingerprint.size_bytes,
                fingerprint.sha256,
                fingerprint.hash_mode,
            )
            for fingerprint in fingerprints
        ]
        return cls(
            name=str(name),
            kind="file" if resolved.is_file() else "directory",
            path=str(resolved),
            content_sha256=sha256_json(entries),
            files=fingerprints,
        )

    @classmethod
    def from_stage_manifest(
        cls, name: str, path: Path | str
    ) -> "InputArtifact":
        manifest = StageManifest.read_json(path)
        return cls(
            name=str(name),
            kind="stage_manifest",
            path=str(Path(path).expanduser().resolve()),
            content_sha256=manifest.content_id,
            files=(),
        )

    @classmethod
    def from_qupath_manifest(
        cls, name: str, path: Path | str
    ) -> "InputArtifact":
        manifest = QuPathImageManifest.read_json(path)
        return cls(
            name=str(name),
            kind="qupath_manifest",
            path=str(Path(path).expanduser().resolve()),
            content_sha256=manifest.content_id,
            files=(),
        )

    @classmethod
    def from_qupath_cohort_manifest(
        cls, name: str, path: Path | str
    ) -> "InputArtifact":
        manifest = QuPathCohortManifest.read_json(path)
        return cls(
            name=str(name),
            kind="qupath_cohort_manifest",
            path=str(Path(path).expanduser().resolve()),
            content_sha256=manifest.content_id,
            files=(),
        )

    def validate_current(self, *, mode: str = "fast") -> None:
        path = Path(self.path)
        if self.kind == "stage_manifest":
            manifest = StageManifest.read_json(path)
            current = manifest.content_id
            if current != self.content_sha256:
                raise StaleArtifactError(
                    f"input {self.name!r} changed "
                    f"(manifest={self.content_sha256[:12]}, current={current[:12]})"
                )
            manifest.validate_current(validation_mode=mode)
            return
        elif self.kind == "qupath_manifest":
            manifest = QuPathImageManifest.read_json(path)
            if manifest.content_id != self.content_sha256:
                raise StaleArtifactError(f"input {self.name!r} manifest changed")
            manifest.validate_current(mode=mode)
            return
        elif self.kind == "qupath_cohort_manifest":
            manifest = QuPathCohortManifest.read_json(path)
            if manifest.content_id != self.content_sha256:
                raise StaleArtifactError(f"input {self.name!r} manifest changed")
            manifest.validate_current(mode=mode)
            return
        else:
            if not path.exists():
                raise StaleArtifactError(f"input artifact is missing: {path}")
            actual_paths = (
                {str(path.resolve())}
                if path.is_file()
                else {
                    str(item.resolve())
                    for item in path.rglob("*")
                    if item.is_file()
                }
            )
            recorded_paths = {fingerprint.path for fingerprint in self.files}
            if actual_paths != recorded_paths:
                raise StaleArtifactError(
                    f"input {self.name!r} file inventory changed"
                )
            for fingerprint in self.files:
                try:
                    fingerprint.validate_current(mode=mode)
                except StaleArtifactError as exc:
                    raise StaleArtifactError(
                        f"input {self.name!r} changed: {exc}"
                    ) from exc
            # Fast/stat validation proves the recorded files are unchanged
            # without rebuilding the aggregate content hash.
            return
        if current != self.content_sha256:
            raise StaleArtifactError(
                f"input {self.name!r} changed "
                f"(manifest={self.content_sha256[:12]}, current={current[:12]})"
            )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "InputArtifact":
        return cls(
            name=str(value["name"]),
            kind=str(value["kind"]),
            path=str(value["path"]),
            content_sha256=str(value["content_sha256"]),
            files=tuple(
                FileFingerprint.from_dict(item) for item in value.get("files", ())
            ),
        )


@dataclass(frozen=True)
class StageManifest:
    """Content-derived contract for one completed pipeline stage."""

    manifest_version: int
    stage: str
    method_version: str
    expected_donors: tuple[str, ...]
    completed_donors: tuple[str, ...]
    inputs: tuple[InputArtifact, ...]
    output: DatasetSnapshot
    sidecars: tuple[FileFingerprint, ...]
    config: ConfigSnapshot
    code: CodeMetadata
    created_at_utc: str
    content_id: str

    @classmethod
    def create(
        cls,
        *,
        stage: str,
        method_version: str,
        expected_donors: Sequence[str],
        inputs: Sequence[InputArtifact],
        output_root: Path | str,
        sidecar_paths: Sequence[Path | str] = (),
        config: Any,
        code: CodeMetadata,
    ) -> "StageManifest":
        if not stage.strip() or not method_version.strip():
            raise ContractError("stage and method_version are required")
        expected = tuple(sorted(map(str, expected_donors)))
        if not expected or len(set(expected)) != len(expected):
            raise ContractError("expected_donors must be a non-empty unique set")
        input_tuple = tuple(sorted(inputs, key=lambda item: item.name))
        if len({item.name for item in input_tuple}) != len(input_tuple):
            raise ContractError("input artifact names must be unique")
        for item in input_tuple:
            item.validate_current()
        output = snapshot_parquet_dataset(
            output_root, expected_donors=expected
        )
        sidecars = tuple(
            sorted(
                (
                    FileFingerprint.capture(path)
                    for path in sidecar_paths
                ),
                key=lambda item: item.path,
            )
        )
        if len({item.path for item in sidecars}) != len(sidecars):
            raise ContractError("stage sidecar paths must be unique")
        draft = cls(
            manifest_version=MANIFEST_VERSION,
            stage=stage.strip(),
            method_version=method_version.strip(),
            expected_donors=expected,
            completed_donors=output.donor_ids,
            inputs=input_tuple,
            output=output,
            sidecars=sidecars,
            config=ConfigSnapshot.capture(config),
            code=code,
            created_at_utc=datetime.now(timezone.utc).isoformat(),
            content_id="",
        )
        return replace(draft, content_id=sha256_json(draft._identity_payload()))

    def _identity_payload(self) -> dict[str, Any]:
        value = self.to_dict()
        value.pop("content_id", None)
        return _content_identity(value)

    def to_dict(self) -> dict[str, Any]:
        return _normalise_json(asdict(self))

    def verify_identity(self) -> None:
        if self.manifest_version != MANIFEST_VERSION:
            raise ContractError(f"unsupported stage manifest version {self.manifest_version}")
        if sha256_json(self._identity_payload()) != self.content_id:
            raise ContractError("stage manifest content_id is invalid")
        if self.completed_donors != self.expected_donors:
            raise PartialArtifactError(
                f"{self.stage}: completed donors do not equal expected donors"
            )

    def validate_current(
        self,
        *,
        config: Any | None = None,
        validate_code: bool = True,
        validate_inputs: bool = True,
        validation_mode: str = "fast",
    ) -> None:
        """Fail if any recorded input/config/code/output is partial or stale."""
        self.verify_identity()
        if validate_inputs:
            for item in self.inputs:
                item.validate_current(mode=validation_mode)
        if config is not None and ConfigSnapshot.capture(config) != self.config:
            raise StaleArtifactError(f"{self.stage}: configuration changed")
        if validate_code:
            current_code = self.code.capture_current()
            if current_code.source_sha256 != self.code.source_sha256:
                raise StaleArtifactError(f"{self.stage}: producing code changed")
        try:
            validate_dataset_snapshot_current(
                self.output, mode=validation_mode
            )
        except StaleArtifactError as exc:
            raise StaleArtifactError(f"{self.stage}: {exc}") from exc
        for sidecar in self.sidecars:
            try:
                sidecar.validate_current(mode=validation_mode)
            except StaleArtifactError as exc:
                raise StaleArtifactError(
                    f"{self.stage}: sidecar changed: {exc}"
                ) from exc

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "StageManifest":
        result = cls(
            manifest_version=int(value["manifest_version"]),
            stage=str(value["stage"]),
            method_version=str(value["method_version"]),
            expected_donors=tuple(map(str, value["expected_donors"])),
            completed_donors=tuple(map(str, value["completed_donors"])),
            inputs=tuple(InputArtifact.from_dict(item) for item in value["inputs"]),
            output=DatasetSnapshot.from_dict(value["output"]),
            sidecars=tuple(
                FileFingerprint.from_dict(item)
                for item in value.get("sidecars", ())
            ),
            config=ConfigSnapshot.from_dict(value["config"]),
            code=CodeMetadata.from_dict(value["code"]),
            created_at_utc=str(value["created_at_utc"]),
            content_id=str(value["content_id"]),
        )
        result.verify_identity()
        return result

    @classmethod
    def read_json(cls, path: Path | str) -> "StageManifest":
        with Path(path).open(encoding="utf-8") as handle:
            return cls.from_dict(json.load(handle))


@dataclass(frozen=True)
class StageReference:
    stage: str
    manifest_path: str
    manifest_content_id: str
    output_content_sha256: str

    @classmethod
    def capture(cls, path: Path | str) -> "StageReference":
        resolved = Path(path).expanduser().resolve()
        manifest = StageManifest.read_json(resolved)
        return cls(
            stage=manifest.stage,
            manifest_path=str(resolved),
            manifest_content_id=manifest.content_id,
            output_content_sha256=manifest.output.content_sha256,
        )

    def validate_current(
        self, *, validate_stage: bool = True, validation_mode: str = "fast"
    ) -> StageManifest:
        manifest = StageManifest.read_json(self.manifest_path)
        if (
            manifest.content_id != self.manifest_content_id
            or manifest.output.content_sha256 != self.output_content_sha256
        ):
            raise StaleArtifactError(f"run stage reference changed: {self.stage}")
        if validate_stage:
            manifest.validate_current(validation_mode=validation_mode)
        return manifest

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "StageReference":
        return cls(
            stage=str(value["stage"]),
            manifest_path=str(value["manifest_path"]),
            manifest_content_id=str(value["manifest_content_id"]),
            output_content_sha256=str(value["output_content_sha256"]),
        )


@dataclass(frozen=True)
class RunManifest:
    """Content-derived collection of exact stage manifests for one run."""

    manifest_version: int
    run_name: str
    expected_donors: tuple[str, ...]
    stages: tuple[StageReference, ...]
    config: ConfigSnapshot
    code: CodeMetadata
    created_at_utc: str
    content_id: str

    @classmethod
    def create(
        cls,
        *,
        run_name: str,
        expected_donors: Sequence[str],
        stage_manifest_paths: Sequence[Path | str],
        config: Any,
        code: CodeMetadata,
    ) -> "RunManifest":
        if not run_name.strip():
            raise ContractError("run_name is required")
        expected = tuple(sorted(map(str, expected_donors)))
        references = tuple(StageReference.capture(path) for path in stage_manifest_paths)
        if not references:
            raise ContractError("a run manifest must reference at least one stage")
        if len({reference.stage for reference in references}) != len(references):
            raise ContractError("a run manifest cannot reference a stage twice")
        for reference in references:
            # A run is a collection of current stage artifacts, not merely a
            # collection of JSON files that existed when it was assembled.
            manifest = reference.validate_current(
                validate_stage=True,
                validation_mode="fast",
            )
            if manifest.expected_donors != expected:
                raise PartialArtifactError(
                    f"stage {manifest.stage} donor contract differs from run donor contract"
                )
        draft = cls(
            manifest_version=MANIFEST_VERSION,
            run_name=run_name.strip(),
            expected_donors=expected,
            stages=references,
            config=ConfigSnapshot.capture(config),
            code=code,
            created_at_utc=datetime.now(timezone.utc).isoformat(),
            content_id="",
        )
        return replace(draft, content_id=sha256_json(draft._identity_payload()))

    def _identity_payload(self) -> dict[str, Any]:
        value = self.to_dict()
        value.pop("content_id", None)
        return _content_identity(value)

    def to_dict(self) -> dict[str, Any]:
        return _normalise_json(asdict(self))

    def verify_identity(self) -> None:
        if self.manifest_version != MANIFEST_VERSION:
            raise ContractError(f"unsupported run manifest version {self.manifest_version}")
        if sha256_json(self._identity_payload()) != self.content_id:
            raise ContractError("run manifest content_id is invalid")

    def validate_current(
        self,
        *,
        config: Any | None = None,
        validate_stages: bool = True,
        validation_mode: str = "fast",
    ) -> None:
        self.verify_identity()
        if config is not None and ConfigSnapshot.capture(config) != self.config:
            raise StaleArtifactError(f"{self.run_name}: run configuration changed")
        current_code = self.code.capture_current()
        if current_code.source_sha256 != self.code.source_sha256:
            raise StaleArtifactError(f"{self.run_name}: run code changed")
        for reference in self.stages:
            manifest = reference.validate_current(
                validate_stage=validate_stages,
                validation_mode=validation_mode,
            )
            if manifest.expected_donors != self.expected_donors:
                raise PartialArtifactError(
                    f"{manifest.stage}: donor set differs from run donor set"
                )

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "RunManifest":
        result = cls(
            manifest_version=int(value["manifest_version"]),
            run_name=str(value["run_name"]),
            expected_donors=tuple(map(str, value["expected_donors"])),
            stages=tuple(StageReference.from_dict(item) for item in value["stages"]),
            config=ConfigSnapshot.from_dict(value["config"]),
            code=CodeMetadata.from_dict(value["code"]),
            created_at_utc=str(value["created_at_utc"]),
            content_id=str(value["content_id"]),
        )
        result.verify_identity()
        return result

    @classmethod
    def read_json(cls, path: Path | str) -> "RunManifest":
        with Path(path).open(encoding="utf-8") as handle:
            return cls.from_dict(json.load(handle))


def _atomic_write_json(path: Path, value: Mapping[str, Any]) -> Path:
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    try:
        temporary.write_text(
            json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


__all__ = [
    "ChannelMapping",
    "CodeMetadata",
    "ConfigSnapshot",
    "ContractError",
    "DatasetSnapshot",
    "DonorDatasetSnapshot",
    "DuplicateObjectError",
    "FileFingerprint",
    "GeometryContractError",
    "InputArtifact",
    "PartialArtifactError",
    "PanelChannelMap",
    "QuPathCohortManifest",
    "QuPathImageManifest",
    "RunManifest",
    "SchemaContractError",
    "StageManifest",
    "StaleArtifactError",
    "canonical_json",
    "dataframe_content_sha256",
    "fingerprint_path",
    "geometry_present",
    "hash_identifiers",
    "image_id_to_donor_id",
    "inspect_qupath_geojson",
    "iter_qupath_geojson_features",
    "sha256_file",
    "sha256_sampled_file",
    "sha256_json",
    "snapshot_parquet_dataset",
    "validate_dataset_snapshot_current",
]
