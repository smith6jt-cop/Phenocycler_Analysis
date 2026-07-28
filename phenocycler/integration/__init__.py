"""
phenocycler.integration — PhenoCycler (protein) <-> Xenium (RNA) integration.

Picks up where the core pipeline stops. ``phenocycler.pipeline`` ends at
``data/phenotype/broad/`` (eight broad lineages, zero Unassigned); this subpackage takes
that output, brings the Xenium side to the same table shape, and integrates them.

Two modes, and the distinction is the whole point:

``sequential``
    Serial sections cut from one FFPE block. The two modalities do **not** share cells, so
    cell-to-cell pairing is invalid no matter how good the registration is. Integration
    happens at four honest levels — donor, islet/structure, niche, and (explicitly labelled
    as inference) pseudo-cell.

``same_slide``
    One physical section imaged twice. Cells *are* the same cells, so a genuine paired
    protein+RNA matrix is recoverable.

``sameslide.py`` and the sequential-only steps refuse to run in the wrong mode. That guard
is a hard error rather than a warning because the validity boundary of every downstream
claim rests on it.

Stage order (``pipeline.py``)::

    manifest      S0   donor <-> Xenium-run pairing, donor_id/roi stamping
    export_pheno  S1a  phenotype/broad + cells + restore_gated  -> cells_pheno/
    import_xenium S1b  h5ad | zarr | bundle                     -> cells_xen/
    vocab         S2   harmonised lineage vocabulary + protein<->gene crosswalk
    structures    S3   tissue / islets / ducts / vessels, same algorithm both sides
    register      S4   coarse rigid -> affine -> non-rigid (image + point-set tracks)
    transform     S5   warp coordinates and masks into the fixed frame
    match         S6   islet <-> islet assignment            [sequential]
    grid          S7   kNN-composition niches + registered-frame hex grid
    donor         S8   registration-free donor-level concordance
    crossmodal    S9   pseudo-cell linking                    [sequential, inference]
    sameslide     S11  cell <-> cell assignment               [same_slide]
    qc            S10  gates, nulls, report

Everything after S1 speaks the modality-agnostic cell-table contract in ``contract.py``,
which is what lets registration, matching and QC be written once for both modalities and
both modes.

The heavy readers (anndata / spatialdata / SimpleITK / opencv) are optional: install with
``pip install -e '.[integration]'`` or use ``environment-integration.yml``. Importing this
package does not pull them in — each module raises a clear error naming the extra only when
the dependency is actually needed, so the core pipeline keeps working in the lean env.
"""

from __future__ import annotations

from .contract import (
    CELL_TABLE_COLUMNS,
    MODALITIES,
    CellTable,
    feature_columns,
    read_cell_table,
    validate_cell_table,
    write_cell_table,
)

__all__ = [
    "CELL_TABLE_COLUMNS",
    "MODALITIES",
    "CellTable",
    "feature_columns",
    "read_cell_table",
    "validate_cell_table",
    "write_cell_table",
]
