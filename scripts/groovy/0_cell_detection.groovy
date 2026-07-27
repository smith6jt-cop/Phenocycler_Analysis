/*
 * 0_cell_detection.groovy
 *
 * Workflow step 1: InstanSeg cell detection over the Tissue annotation ONLY.
 *
 * WHY THIS SCRIPT EXISTS
 * ----------------------
 * Detection used to be a manual GUI step — the only unscripted stage in the pipeline — and it
 * silently corrupted 8 of 20 donors. The QuPath workflow history shows the affected runs used:
 *
 *     selectAnnotations()            <-- selects Tissue AND every nested Islet_N child
 *     InstanSeg...detectObjects()
 *
 * InstanSeg's detectObjects() runs inference once per SELECTED parent object and does not
 * deduplicate across parents (PruneObjectOutputHandler dedupes only across tiles WITHIN one
 * parent). So every nucleus inside an islet was segmented twice — once under Tissue, once under
 * its islet. The two passes tile from different ROI origins, so the copies are near-identical but
 * not equal: centroids ~1 px apart, areas within a fraction of a percent, different vertex counts,
 * two distinct UUIDs. 8_hierarchy_parenting.groovy then re-parents both to the islet, which is why
 * 100% of the duplicates landed in islet cores.
 *
 * Measured damage: 155,669 duplicate pairs across donors 6450, 6505, 6523, 6534, 6566, 6591, 6623.
 * Roughly half to two-thirds of every islet cell in those donors was a duplicate. Nothing
 * downstream detected it — the UUIDs are distinct, so every uniqueness check passed.
 *
 * THE RULES THIS SCRIPT ENFORCES
 * ------------------------------
 *   1. NEVER selectAnnotations(). Select the Tissue annotation explicitly, by classification.
 *   2. Refuse to run if any non-Tissue annotation is selected, so the failure cannot recur.
 *   3. clearDetections() first, so a re-run replaces rather than accumulates.
 *   4. Verify afterwards that no two detections share a centroid — fail loudly if they do.
 *
 * Run it BEFORE 1_panc_objects.groovy (which creates the Islet annotations). The rule-2 guard makes
 * it safe to run afterwards too, but running first keeps the ordering unambiguous.
 *
 * Headless:
 *   /home/smith6jt/QuPath/bin/QuPath script \
 *       --project /home/smith6jt/IO60panc2nd/project.qpproj \
 *       --image "6539_Scan1.er.qptiff - resolution #1" \
 *       /home/smith6jt/IO60panc2nd/scripts/0_cell_detection.groovy
 */

import qupath.lib.images.servers.ColorTransforms

// ==================== PARAMETERS ====================
// Matches every run recorded in the project workflow history (QuPath 0.7.0).
String MODEL_PATH   = "/home/smith6jt/QuPath/v0.7/instanseg/downloaded/fluorescence_nuclei_and_cells-0.1.1"
String DEVICE       = "gpu0"      // "cpu" if no CUDA device
int    TILE_DIMS    = 512
int    PAD          = 32
int    N_THREADS    = 32
boolean MAKE_MEASUREMENTS = false // measurements come from 2_fast_cell_measurements.groovy
String TISSUE_CLASS = "Tissue"
double DUP_RADIUS_UM = 1.0        // post-run duplicate check; matches phenocycler.cell_qc
double DUP_AREA_TOL  = 0.15
// ====================================================

def server = getCurrentServer()
if (server == null) { println "ERROR: no image open."; return }
def imageName = getProjectEntry()?.getImageName() ?: server.getMetadata()?.getName()
println "[detect] image: ${imageName}"

// ---- 1. select ONLY the tissue annotation ---------------------------------
resetSelection()
selectObjectsByClassification(TISSUE_CLASS)
def selected = getSelectedObjects()
if (selected.isEmpty()) {
    println "ERROR: no '${TISSUE_CLASS}' annotation found. Create it before running detection."
    return
}

// ---- 2. refuse anything that is not Tissue --------------------------------
// This is the guard for the exact failure above: if an Islet (or any other) annotation is in the
// selection, InstanSeg would segment its cells a second time.
def offending = selected.findAll { it.getPathClass()?.toString() != TISSUE_CLASS }
if (!offending.isEmpty()) {
    println "ERROR: selection contains ${offending.size()} non-${TISSUE_CLASS} object(s):"
    offending.groupBy { it.getPathClass()?.toString() ?: 'Unclassified' }
             .each { k, v -> println "         ${k}: ${v.size()}" }
    println "       Detection would segment their contents a SECOND time. Refusing to run."
    return
}
println "[detect] selected ${selected.size()} '${TISSUE_CLASS}' annotation(s) — nothing else."

// ---- 3. clear existing detections so a re-run replaces, never accumulates --
int before = getDetectionObjects().size()
if (before > 0) {
    println "[detect] clearing ${before} existing detections"
    clearDetections()
    resetSelection()
    selectObjectsByClassification(TISSUE_CLASS)
}

// ---- 4. run InstanSeg on every channel present ----------------------------
def presentNames = server.getMetadata().getChannels()
    .collect { it?.getName() }.findAll { it && it.trim() }
def extractors = presentNames.collect { ColorTransforms.createChannelExtractor(it) }
println "[detect] ${presentNames.size()} input channels, model ${MODEL_PATH.split('/')[-1]}"

long t0 = System.currentTimeMillis()
qupath.ext.instanseg.core.InstanSeg.builder()
    .modelPath(MODEL_PATH)
    .device(DEVICE)
    .inputChannels(extractors)
    .tileDims(TILE_DIMS)
    .interTilePadding(PAD)
    .nThreads(N_THREADS)
    .makeMeasurements(MAKE_MEASUREMENTS)
    .randomColors(false)
    .outputType("Default")
    .build()
    .detectObjects()
fireHierarchyUpdate()

def cells = getDetectionObjects()
println String.format("[detect] %,d detections in %.0fs", cells.size(),
        (System.currentTimeMillis() - t0) / 1000.0)

// ---- 5. verify no duplicates ----------------------------------------------
// Grid-bucket the centroids so this is O(n) rather than O(n^2) on ~2M cells.
def cal = server.getPixelCalibration()
double pxSize = cal.hasPixelSizeMicrons() ? cal.getAveragedPixelSizeMicrons() : 1.0
double rPx = DUP_RADIUS_UM / pxSize
def buckets = [:].withDefault { [] }
for (c in cells) {
    def roi = c.getROI()
    if (roi == null) continue
    buckets[[(int) (roi.getCentroidX() / rPx), (int) (roi.getCentroidY() / rPx)]] << c
}
int dupPairs = 0
for (entry in buckets) {
    def (bx, by) = entry.key
    def near = []
    for (dx in -1..1) for (dy in -1..1) near.addAll(buckets[[bx + dx, by + dy]])
    for (a in entry.value) {
        def ra = a.getROI()
        for (b in near) {
            if (System.identityHashCode(b) <= System.identityHashCode(a)) continue
            def rb = b.getROI()
            double d = Math.hypot(ra.getCentroidX() - rb.getCentroidX(),
                                  ra.getCentroidY() - rb.getCentroidY())
            if (d > rPx) continue
            double aa = ra.getArea(), ab = rb.getArea()
            if (Math.abs(aa - ab) / Math.max(Math.max(aa, ab), 1e-9) < DUP_AREA_TOL) dupPairs++
        }
    }
}
if (dupPairs > 0) {
    println "ERROR: ${dupPairs} duplicate detection pair(s) found after detection."
    println "       Two objects share a centroid within ${DUP_RADIUS_UM} um with matching area."
    println "       This is the selectAnnotations() failure — the detections are NOT usable."
} else {
    println "[detect] duplicate check: PASS (0 pairs within ${DUP_RADIUS_UM} um)"
}
println "[detect] Done."
