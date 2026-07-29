# Running this on HiPerGator — setup and the basic analyses

A working outline for the 20-donor cohort: PhenoCycler protein and Xenium 5K RNA, pancreas
and pancreatic lymph node, on UF Research Computing.

It is written for the state you are actually in — refining REDSEA→RESTORE marker positivity —
so [Phase 4](#phase-4--the-positivity-refinement-loop) is the centre of gravity, not an
afterthought. Everything before it exists to get you to a cohort you can iterate on;
everything after it depends on positivity being settled, and quietly goes stale when it
changes.

Commands are the real ones (checked against `argparse`); paths are the real ones (checked
against `PipelineConfig`). Where a number is yours to supply rather than the repo's, it says so.

---

## The shape of it

```
OFF-CLUSTER                          ON HIPERGATOR
────────────                         ─────────────
QuPath segmentation                  [1] cells_parquet    DuckDB, one job
  ├─ Cellmeasurements.csv  ────────▶ [2] redsea           SLURM array, one SECTION per task
  └─ cells__<image>.geojson ───────▶ [3] restore + extra  one job, cohort-wide fit
                                     [4] hormone_floor ─┐
Xenium runs                          [5] lineage        ├─ one job
  └─ output bundles ───────┐        [6] qupath          │
                            │        [7] figures        ┘
                            │        [8] reassess_diag  ← the acceptance test
                            │              │
                            │              └──▶ iterate on positivity ──┐
                            │                                           │
                            └──────────────▶ [9] integration ◀──────────┘
                                              (second conda env)
```

Two units matter and are easy to confuse:

- **Section** — one image, one GeoJSON, one RESTORE threshold set, one coordinate frame.
  `6539` (pancreas) and `6539pln` (lymph node) are *two* sections. This is the pipeline's
  work unit and what sizes the SLURM array.
- **Donor** — `6539`. This is what joins to metadata and to the Xenium manifest.

`cfg.discover_sections()` and `cfg.discover_donors()` are both real methods and return
different lists. 20 donors × 2 tissues is ~40 sections.

---

## Phase 0 — what happens off-cluster

Neither of these is in this repo, and neither can be done from a SLURM job.

**0.1 QuPath, per image.** Segment, then export two things:
- the per-cell measurement table → one `Cellmeasurements.csv` for the cohort;
- full-resolution boundaries → `cells__<image>.geojson`, one per image, via
  `scripts/groovy/export_cells_geojson.groovy`.

Edit the script's `folder`/`outDir` first — it hardcodes `~/Phenocycler_Analysis/data/...`,
which is not where your data will live on HiPerGator. The filename matters: REDSEA rebuilds
it from the section's `Image` value with `[^A-Za-z0-9._-]` → `_`, so transfer the files with
names byte-preserved.

**0.2 Xenium bundles.** You need the output bundle per run (`cells.parquet`,
`cell_boundaries.parquet`, transcripts) — the polygons ship there and *not* in the zarr or a
slimmed h5ad. If you have phenotyped h5ads from `Xenium_Analysis`, note where their
`data/processed` root is; Phase 5 needs it.

**0.3 Decide the section naming for the post-Xenium re-stains — before segmentation.**
Those qptiffs go through the core pipeline like any other image, so their filename has to
carry a token `phenocycler/sections.py` recognises: `xen`, `xpanc`, `xenpln`, `xpln`
(e.g. `6539xen_Scan1.er.qptiff`). An unrecognised token is **refused**, not guessed at. Renaming
after segmentation means re-exporting every GeoJSON, so agree the convention now.

---

## Phase 1 — one-time setup on HiPerGator

**1.1 Clone with submodules.**

```bash
git clone --recurse-submodules https://github.com/smith6jt-cop/Phenocycler_Analysis.git
cd Phenocycler_Analysis
git submodule update --init --recursive
ls external/RESTORE/python_code            # must exist, or step 3 dies in seconds
ls external/XeniumPanelExplorer/data/panel_roles.csv
```

Both `external/` checkouts are load-bearing and both fail late and confusingly if missing:
RESTORE is the normalization itself; XeniumPanelExplorer supplies the panel taxonomy the
vocabulary crosswalk is derived from. Verify them now, not after a 24-hour job.

**1.2 Environments.** Two, deliberately — the core pipeline stays lean so it does not inherit
the Xenium stack.

```bash
bash setup.sh                                  # submodules + phenocycler_analysis env + verify
conda activate phenocycler_analysis
pip install spams-bin                          # SSC model for RESTORE, if not already pulled

# Only needed for Phase 5:
conda env create -f environment-integration.yml
conda activate phenocycler_integration && pip install -e .
```

`setup.sh` does *not* create the integration env — that line is yours. Note also that conda
envs default to `~/.conda/envs`, which is on your home quota; `conda config --add envs_dirs
/blue/<group>/<user>/envs` first if that is tight.

**1.3 Put the data somewhere with room.** `config.ini` ships `data_dir = data`, i.e. inside
the repo, which on HiPerGator means your home quota. Redirect it:

```bash
export PHENOCYCLER_DATA_DIR=/blue/<group>/<user>/phenocycler/data
export PHENOCYCLER_IMAGES_DIR=/blue/<group>/<user>/phenocycler/images
export PHENOCYCLER_CELLS_CSV=/blue/<group>/<user>/phenocycler/Cellmeasurements.csv
export PHENOCYCLER_DONOR_METADATA=/blue/<group>/<user>/phenocycler/donor_metadata_panc.xlsx
export TMPDIR=/blue/<group>/<user>/tmp            # DuckDB spills here in step 1
```

Or edit `[paths]`. Prefer the environment variables while you are iterating — see §4.6 for
why editing `config.ini` has a side effect.

Budget: the QuPath CSV is ~138 GB; the parquet cohort is tens of GB; and `--keep-mask
--save-intermediates` (on by default in the array, see §4.3) adds a full-resolution int32 mask
plus an `.npz` per section on top of that — plausibly hundreds of GB across ~40 sections. Drop
those two flags if you are not sweeping REDSEA parameters.

**1.4 Fill in the SBATCH placeholders.** Three literals appear in all five job scripts:
`your_email@ufl.edu`, `your_qos`, `your_account`. For account and QOS you can skip the editing
and `export SBATCH_ACCOUNT=... SBATCH_QOS=...` instead — sbatch input environment variables
take precedence over in-script `#SBATCH` directives. There is no equivalent for `--mail-user`,
so that one still needs editing or deleting.

**1.5 Submit from the repo root, always.** Every `#SBATCH --output` is the relative path
`logs/...`, and SLURM opens those files *before* the job's own `mkdir -p logs` runs. Submitting
from elsewhere kills the job before Python starts, with nothing in any log to say why.
`mkdir -p logs` once, up front.

---

## Phase 2 — smoke test before you commit to the queue

All of this runs on a login or dev node in minutes. It is the only way to find out that a
submodule, `spams`, or a config path is wrong before a 12-hour array discovers it for you.

```bash
pytest tests/ -q                                     # ~30 s, fully synthetic, no data needed
python -m phenocycler.pipeline --status              # renders with no data; proves config loads
python -m phenocycler.integration.pipeline --status
python -m phenocycler.integration.manifest --report  # donor <-> Xenium pairing table

python -c "import spams; print('spams ok')"
python -c "from phenocycler import load_config; c=load_config(); print(c.data_dir, c.cells_csv)"
```

Then one real section, end to end, timed — this is what sizes the array:

```bash
python -c "
from phenocycler import load_config, cells_parquet
cells_parquet.build_cells_parquet(load_config(), limit=200_000)"

/usr/bin/time -v python -m phenocycler.redsea --donor <section_key> 2>&1 | tee logs/smoke.log
```

Read the wall time and peak RSS out of that log and size `--time` and `--mem` from them. The
shipped `12:00:00` / `64gb` are guesses, and REDSEA's peak is driven by your image dimensions:
the full-resolution int32 mask, plus `distance_transform_edt` returning a float64 distance
array *and* a 2×int64 index array over the same `(H, W)`. Print one section's shape before
believing 64 GB.

---

## Phase 3 — the protein pipeline

```bash
bash scripts/slurm/run_full_pipeline.sh 40     # 40 SECTIONS -> array 0-39
```

That submits `cells → redsea[array] → restore(+extra) → floor/lineage/qupath/figures` as an
`afterok` chain. Stage by stage:

| # | stage | writes | notes |
|---|---|---|---|
| 1 | `cells` | `data/cells/donor_id=*` | DuckDB, out-of-core; the only stage that reads the CSV |
| 2 | `redsea` | `data/cells_redsea/donor_id=*/data_0.parquet` | one array task per section |
| 3 | `restore` | `data/restore_gated_redsea/`, `data/restore_thresholds_redsea.csv` | cohort-wide fit, per-donor apply |
| 4 | `restore_extra` | `data/restore_gated_redsea_extra/` | CD99 / B3TUBB / MPO, separate so the validated 10-marker gates stay byte-identical |
| 5 | `hormone_floor` | `data/restore_gated_redsea/` (in place) | backs up to `data/restore_gated_redsea.pre_hormonefloor/` first |
| 6 | `lineage` | `data/phenotype/broad/donor_id=*` | 7 classes, zero Unassigned |
| 7 | `qupath` | `data/phenotype/qupath_class/pheno_class_*.csv` | UUID-keyed, for re-import |
| 8 | `figures` | `data/phenotype/celltype_marker_{dotplot,heatmap}.png` | identity QC |

### 3.1 Do not trust `afterok` alone

`afterok` checks the job's exit code, and per-section failures inside REDSEA and RESTORE are
**logged, not raised** (`parallel.map_donors` defaults to `on_error="log"`). Separately,
`pipeline --status` reports a stage `OK` if it wrote *at least one* partition. So one surviving
section makes a stage look complete, and every later run skips it.

Count partitions between stages:

```bash
for d in cells cells_redsea restore_gated_redsea restore_gated_redsea_extra; do
  printf "%-28s %s\n" "$d" "$(ls -d "$PHENOCYCLER_DATA_DIR/$d"/donor_id=* 2>/dev/null | wc -l)"
done
grep -c FAILED logs/redsea_*.out
```

All four numbers should equal your section count. This is the single most valuable check in
the whole document.

### 3.2 What the pipeline actually decides

Broad lineage is a strict hierarchy, not a classifier — `phenocycler/lineage.py`:

```
Endocrine     INS|GCG|SST positive, or bright CD99          (highest precedence)
Immune        CD3e|CD20|CD163|MPO positive                   MPO = neutrophils, folded in
Endothelial   CD31 positive
Neural        B3TUBB positive
─ structural background: argmax over Pan_Cytokeratin / Vimentin / SMA _norm
Epithelial    everything below all three structural gates    → flagged `epi_default`
```

Seven classes, no Unassigned bucket. `epi_default` marks cells that landed in Epithelial by
exhaustion rather than by evidence — it is the honest sink, and worth tracking as a fraction.

**There are no fine/specific cell types today.** `endocrine_subtype` (Beta/Alpha/Delta from
which hormone is brightest) is derived in the *integration* layer, not the core pipeline, and
nothing else subtypes. When you get to specific types, that is new work, not a flag to flip.

---

## Phase 4 — the positivity refinement loop

This is where you are. The chain is:

```
REDSEA           alpha, comp_mode, edge_radius, gap_bridge    corrects raw MFI
   ↓
RESTORE          per-image threshold per marker (SSC)         writes _pos, _norm, _log2r
   ↓                                                          _norm = raw / threshold
hormone_floor    _pos := _norm >= K, for 6 markers            K = 5 hormone, 2 immune
   ↓
lineage          hierarchy over _pos                          + CD99 and MPO gated in-script
```

`_norm` is threshold-relative by construction, which is what makes a floor meaningful: `K=5`
means "five times the per-image RESTORE threshold", not an absolute intensity.

### 4.1 Which knob lives where — this is the part that wastes days

| knob | applied in | re-run to change it |
|---|---|---|
| `alpha`, `comp_mode`, `edge_radius`, `gap_bridge` | `redsea` | everything downstream |
| `model`, `subsample`, `robust_factor`, marker pairs | `restore` | restore → floor → lineage |
| `hormone_min_norm` (INS/GCG/SST) | `hormone_floor` | floor → lineage |
| `immune_min_norm` — CD3e/CD20/CD163 | `hormone_floor` | floor → lineage |
| `immune_min_norm` — **MPO** | **`lineage`** | lineage only |
| `cd99_bright` | **`lineage`** | lineage only |
| `B3TUBB` positivity | RESTORE only, no floor anywhere | restore --extra |

Two of those are traps. Changing `immune_min_norm` and re-running only `hormone_floor` leaves
**MPO** gated at the old K, because MPO lives in the extra directory and is floored inside
`lineage` instead. And changing `cd99_bright` and re-running `hormone_floor` does **nothing at
all** — it is only ever read by `lineage`.

### 4.2 The cheapest correct re-run

For a K change, on a dev node, two or three donors:

```bash
python -m phenocycler.hormone_floor \
    --gated-dir "$PHENOCYCLER_DATA_DIR/restore_gated_redsea.pre_hormonefloor" \
    --out-dir   "$PHENOCYCLER_DATA_DIR/restore_gated_redsea" \
    --donors 6539 6476 6579 --jobs 8 | tee logs/floor_K.log

python -m phenocycler.lineage --donors 6539 6476 6579 --jobs 8
python -m phenocycler.reassess_diag --donors 6539 6476 6579
```

Minutes, not hours. Flooring **from the pre-floor backup** is the important part: flooring an
already-floored table is not wrong (the floor recomputes `_pos` from `_norm`, which it never
touches) but starting from clean RESTORE gates makes the result reproducible from a known state.

Note the module CLIs take `--donors`; `python -m phenocycler.pipeline` does **not**. Per-donor
scoping only exists at module level.

Through the orchestrator, cohort-wide:

```bash
python -m phenocycler.pipeline --only hormone_floor lineage qupath figures --force
```

**`--force` is not optional here.** `lineage` is the one stage with a custom staleness check,
and it only asks whether the columns have the current 7-class shape — not what K they were
computed at. Without `--force` you get freshly floored gates and the *old* lineage calls, and
the status table says `OK (7-class)`. This is the most likely way a threshold change appears
to do nothing.

### 4.3 Sweeping REDSEA without re-reading qptiffs

The array passes `--keep-mask --save-intermediates`, which persists the rasterized mask and
the `data`/`edge`/`sizes`/`contact` arrays per section. That makes a parameter sweep cheap:

```bash
python -m phenocycler.redsea --donor <section> --from-intermediates --alpha 0.8
python -m phenocycler.redsea --donor <section> --from-intermediates --comp-mode 1
python -m phenocycler.redsea --donor <section> --from-intermediates \
    --alpha-per-channel 'Pan_Cytokeratin=2,Vimentin=1.5'
```

Seconds per section instead of a full qptiff read and re-rasterization. This is the entire
reason those two flags are on by default; if you are not sweeping, drop them and save the disk.

### 4.4 Three rules about re-running RESTORE

1. **Delete the pre-floor backup whenever RESTORE re-runs.** `pipeline`'s floor stage snapshots
   `restore_gated_redsea/` → `restore_gated_redsea.pre_hormonefloor/` only if the backup is
   absent, and then always floors *from* it. After a new RESTORE fit that backup is stale, and
   the next floor will quietly overwrite your new gates with the old ones.
   `rm -rf "$PHENOCYCLER_DATA_DIR/restore_gated_redsea.pre_hormonefloor"` — and **only** then.
   Never delete it for a K-only change; that is the whole point of it.

2. **Never run `--only restore` after a floor.** `restore` and `hormone_floor` share an output
   directory, so re-running RESTORE reverts the six floored `_pos` columns to `_norm >= 1` and
   nothing detects it. Re-run the floor immediately after.

3. **A subset fit is not a subset of the full fit.** RESTORE's robust guard imputes degenerate
   thresholds from the cohort median *of whatever was fitted in that run*, so
   `restore --donors 6539 6476` gives those donors different thresholds than the full cohort
   would. Use `--donors` for inspection; validate a K against a full-cohort fit. `--reuse-threshs`
   (re-apply existing thresholds) and `--skip-apply` (fit and QC only) are the honest cheap paths.

### 4.5 Know when you are done

`phenocycler/reassess_diag.py` calls itself the acceptance yardstick for any change to REDSEA,
and it is the only thing in the repo that answers "did that help?". It is now wired into
`04_lineage.sh`; run it by hand during iteration.

```bash
python -m phenocycler.reassess_diag --donors 6539 6476 6579
# -> data/redsea_reassess/{false_endocrine,realvsfalse,guardrails,endocrine_by_stratum}.csv
```

Objective: endocrine keratin/Vimentin load and impossible co-expression should fall.
Guardrails, which must **not** move: acinar keratin retention, vessel Vimentin, the ND→T1D
β-loss and T-cell-gain trends. A change that improves the objective while breaking a guardrail
has not fixed spillover, it has deleted real signal.

`hormone_floor` itself has no persisted output — its whole verification surface is a stdout
line per donor (`[6539] 412,338 cells | INS>=5 51,204->6,881 | ...`), and `--status` shows it
`OK` merely because the gated directory exists. `tee` the log, or spot-check a partition:

```python
import pandas as pd
df = pd.read_parquet(f"{DATA}/restore_gated_redsea/donor_id=6539/data_0.parquet")
assert (df.INS_pos == (df.INS_norm >= 5)).all()
```

### 4.6 Two things to know before you edit `config.ini`

- `tests/test_config.py::test_scientific_defaults` loads the repo's `config.ini` and pins nine
  values: `redsea_comp_mode`, `redsea_alpha`, `redsea_edge_radius`, `restore_model`,
  `restore_robust`, `restore_robust_factor`, `restore_min_cell_area`, **`hormone_min_norm == 5.0`**
  and **`cd99_bright == 3.0`**. Editing any of those in the file turns the test suite red —
  which is the point, but it means a sweep should go through
  `PHENOCYCLER_HORMONE_MIN_NORM` / `PHENOCYCLER_IMMUNE_MIN_NORM` and the file should change once,
  with the test, when you settle. `immune_min_norm` is *not* pinned, so it can be edited freely.
  There is **no** environment override for `cd99_bright` or for any `[redsea]`/`[restore]` key.
- `cd99_bright` is defined in three places — `config.cd99_bright` drives the lineage call,
  `marker_taxonomy.CD99_BRIGHT` the dotplot, `reassess_diag.CD99_BRIGHT` the diagnostic.
  Changing the ini alone makes the three disagree; a test guards the two module constants.

### 4.7 Freeze the cohort before you settle K

RESTORE fits per-image thresholds and imputes outliers against a cohort median. Adding the
not-yet-exported pLN and post-Xenium sections to `data/cells` and re-running RESTORE therefore
re-fits that median and can shift `_norm` for **every existing section**. A K validated today
is not automatically the right K after those sections land. Either freeze the section list
before tuning, or plan to re-validate after the import.

### 4.8 One inconsistency in the repo worth resolving

Donor **6476** is described two ways. `config.py` and `hormone_floor.py` treat its ~40% CD3e as
the threshold over-call that motivates the `K=2` immune floor. `reassess_diag.py` treats 6476
as the ND-pancreatitis donor whose extreme CD3e is *real* and excludes it from the ND→T1D trend.
Both cannot be right, and which one is decides whether that donor's T cells survive the floor.
You are the one who can settle it; nothing in the code can.

---

## Phase 5 — integrating with Xenium

Second environment, second orchestrator, same idempotent shape.

**5.1 Build the manifest with the processed root — do this first, by hand.**

```bash
conda activate phenocycler_integration
python -m phenocycler.integration.manifest --report          # look before writing
python -m phenocycler.integration.manifest --build \
       --xenium-processed /path/to/Xenium_Analysis/data/processed
```

This matters more than it looks. The `manifest` *stage* inside the pipeline runs without
`--xenium-processed`, so `xenium_h5ad`/`xenium_zarr` stay blank and `import_xenium` falls back
to the raw bundle — which carries coordinates and morphology but **no cell-type labels**. The
run completes, the crosswalk is empty, and donor concordance compares nothing. Build the
manifest yourself first; the stage will then skip.

Two data-entry prerequisites the report will show you:

- `data/integration/donor_overrides.csv` has two rows (`0041323__panc`, `0041326__panc`) with a
  **blank** `donor_id`. They stay `donor_unknown` until a human fills them in.
- `data/integration/xenium_paths.csv` carries absolute `/orange/brusko/xenium/...` paths and
  lists 26 donors, of which your cohort is 20 — and which 20 is recorded nowhere. Verify one
  bundle resolves (`ls <path>/cell_boundaries.parquet`) or set `[integration] xenium_root` to
  rewrite the prefix.

**5.2 Run it.**

```bash
sbatch scripts/slurm/06_integration.sh                 # whole cohort, both tissues
ROI=pln_2 sbatch scripts/slurm/06_integration.sh       # one ROI, to isolate a problem
```

Stages, in sequential mode: `manifest → export_pheno → import_xenium → postxen → vocab →
structures → register → match → cellmatch → grid → donor → crossmodal → qc → figures`.

**5.3 Read the report, not the exit code.** `data/integration/qc/integration_report.txt` grades
each donor's registration:

- **PASS** — both gates met: tissue Dice ≥ 0.80 and islet RMSE ≤ 150 µm.
- **WARN** — exactly one gate met.
- **FAIL** — neither. Excluded from structure, niche and pseudo-cell analyses, but still
  contributes to the registration-free donor-level comparison.

For pancreatic lymph node, islet RMSE is dropped from the decision entirely rather than counted
as a failure — a lymph node has no islets, so there is no landmark to measure against, and it
is graded on Dice alone.

**5.4 Expect `n/a` for pLN on three stages.** `structures`, `match` and `cellmatch` need a
structural unit, and `STRUCTURE_KINDS["pancreatic_lymph_node"]` is empty by design. The line
`[n/a] ... no selected tissue has a structural unit` is the correct output, not a failure.

**5.5 The floor change propagates here too, and the orchestrator will not notice.**
`export_pheno` re-derives `endocrine_subtype` from `cfg.hormone_min_norm` and re-applies
`cfg.cd99_bright`; `structures` seeds islet DBSCAN at `hormone_min_norm`. But every stage skips
when its output glob already matches. After any positivity change:

```bash
python -m phenocycler.integration.pipeline --force
```

Otherwise islets, niches, donor composition and crossmodal anchors keep the old thresholds
indefinitely. `data/cells/` and `data/cells_redsea/` are the only things invariant to every
positivity knob — which is exactly why a blanket `--force` on the *core* pipeline is the wrong
reflex: it re-parses the 138 GB CSV for a threshold experiment.

---

## Phase 6 — what cannot run yet, and why that is fine

| | status | unblocked by |
|---|---|---|
| `postxen` (S1c) | no inputs | post-Xenium qptiffs through QuPath → the full core pipeline, under a recognised section token (§0.3) |
| `xenium_hormone_if` | must stay `false` | `postxen` actually stamping `endocrine_subtype`. Setting it early asserts a hormone-protein measurement you do not have |
| `cellmatch` (S6b) | runs, provisional | today the Xenium side's Beta/Alpha/Delta comes from surrogate cores (PDX1/ISL1/NEUROD1/ABCC8) because INS/GCG/SST are off the 5K panel. Label its output provisional until `postxen` lands |
| `sameslide` (S11) | not applicable | nothing — this cohort is serial-section. The genuine same-slide relationship is `postxen`, which runs in *both* modes on purpose. Do not flip `mode = same_slide` to get it; that breaks `match`/`cellmatch`/`crossmodal` and asserts cell identity across serial sections |

One inference to check: a bare `pLN` section token maps to `pln_1`. That mapping was reasoned
from the 26-donor `xenium_paths.csv`, in which donors 6457, 6534 and 6579 have more than one
pLN Xenium run. If any of those donors' single PhenoCycler lymph-node section is really their
`pln_2`, registration will simply not converge. The fix is a per-donor entry in
`sections.REGION_TO_ROI`, not a looser mapping.

---

## Appendix A — decide vs. leave alone

**Leave alone** — `config.ini` says these "should not usually change", and the starred ones are
additionally pinned by `tests/test_config.py`: `[redsea] alpha=1.0`\*, `comp_mode=0`\*,
`edge_radius=0`\*, `gap_bridge=1`; `[restore] model=SSC`\*, `robust=true`\*,
`robust_factor=3.0`\*, `min_cell_area=5.0`\*, `subsample=15000`, `seed=0`; the marker pairs;
`LINEAGES`; and the integration QC gates (`qc_tissue_dice_min=0.80`, `qc_islet_rmse_max_um=150`).

**Yours to decide:** `data_dir` location; `n_jobs`/`duckdb_threads` against `--cpus-per-task`;
array size and walltime; whether to keep masks and intermediates; whether the cohort is frozen
before tuning K; the K values themselves; and `[integration] mode` (leave `sequential`).

**One nobody thinks of:** `restore_extra` does **not** need re-running when K changes. K
operates on `_norm`, which the extra pass does not touch, and `run_restore_extra` hardcodes
`robust=False` so `robust_factor` changes never reach it either.

## Appendix B — silent-failure checklist

Run these between phases. Each corresponds to a failure that produces no error.

```bash
# 1. Every section made it through every stage (§3.1)
for d in cells cells_redsea restore_gated_redsea restore_gated_redsea_extra; do
  printf "%-28s %s\n" "$d" "$(ls -d "$PHENOCYCLER_DATA_DIR/$d"/donor_id=* | wc -l)"; done

# 2. The floor actually applied at the K you set (§4.5)
grep -h "INS>=" logs/lineage_*.out | head

# 3. Lineage is current, not stale-but-shaped-right (§4.2)
python -m phenocycler.pipeline --status      # must read "OK (7-class)", after --force

# 4. Xenium arrived with labels, not coordinates only (§5.1)
grep -i "importing coordinates only\|donor_unknown" logs/integration_*.out

# 5. Registration verdicts, not the job exit code (§5.3)
cat "$PHENOCYCLER_DATA_DIR/integration/qc/integration_report.txt"

# 6. Crossmodal did not silently lose its anchors
grep -i "anchor" logs/integration_*.out
```

`crossmodal` refuses below `crossmodal_min_anchors = 8` verified protein↔gene anchors, and the
usable set is thinner than the 46 nominal 1:1 pairs because INS, GCG, SST and Vimentin are off
the Xenium panel. Lowering that number is a scientific concession, not a workaround.

---

## Appendix C — what this document changed in the repo

Writing it surfaced five things that would have broken the run it describes. All are fixed:

1. `environment.yml` had no `statsmodels`, which `phenocycler/lineage.py` imports at module
   load. An env built from that file could not run step 5 — and `setup.sh`'s own verification
   imports `lineage`, so setup failed at `[3/4]`.
2. The SLURM chain went `restore → lineage`, skipping `restore --extra` and the floor.
   `lineage` raises `FileNotFoundError: extra-gated markers missing` without the extra pass, so
   step 4 died for every donor after step 3 had already spent its 24 hours.
3. `02_redsea_array.sh` indexed `discover_donors()` while REDSEA consumes sections, so every
   `*pln` section was silently never spillover-corrected and the job still exited 0.
4. `06_integration.sh` defaulted `ROI=panc`, and `--roi` overrides `--tissue` — the shipped
   script ran pancreas only and never touched the lymph nodes.
5. All five job scripts called `conda activate` without sourcing conda's shell hook, which
   under `set -e` kills a non-interactive SLURM job before any Python runs.

Also: `xenium_hormone_if` is now documented in `config.ini` (it was read by the loader but
mentioned nowhere), and the stale "eight lineages" claims left over from the Neutrophil→Immune
merge are corrected in `README.md`, `DATA_README.md`, `pyproject.toml` and
`import_broad_lineage.groovy`.
