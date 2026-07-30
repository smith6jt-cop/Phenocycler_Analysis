/*
 * export_cells_geojson.groovy
 *
 * Per-cell geometry export for the manifest-driven Phenocycler workflow.
 * Exports detections rather than annotations and excludes measurement values;
 * the contracted measurement CSV remains their sole source.
 *
 * Each exported feature carries:
 *   - "id"        : the QuPath detection UUID == object_id in data/cells/*.parquet
 *   - "geometry"         : the cell polygon, in FULL-RES pixel coordinates (downsample 1)
 *   - "nucleusGeometry"  : the nucleus polygon, required by geometry QC and marker localization
 *
 * Output (one file per image):
 *   <OUT_DIR>/cells__<imageTag>.geojson
 *
 * Run headless with an explicit output directory:
 *   QuPath script --project project.qpproj --image "exact image name" \
 *     scripts/groovy/export_cells_geojson.groovy --args "/output/geometry"
 */

String outDir = (
    args != null && args.size() > 0 ? args[0]?.trim() : ""
)
if (outDir == null || outDir.isEmpty()) {
    throw new IllegalArgumentException(
        "Pass the geometry output directory as the first script argument."
    )
}

def server = getCurrentServer()
if (server == null) {
    throw new IllegalStateException("No image is open; geometry export cannot run.")
}

// Cells = detection objects with a ROI. (Annotations like Islet/Tissue are skipped.)
def cells = getDetectionObjects().findAll { it.getROI() != null }
if (cells.isEmpty()) {
    throw new IllegalStateException("No detection objects (cells) found in this image.")
}

// REDSEA/QC require a one-to-one cell+nucleus contract. QuPath's GeoJSON exporter writes the
// nucleus as `nucleusGeometry`, but only when the detection is a PathCellObject with a nucleus ROI.
def missingNucleus = cells.findAll {
    !(it.metaClass.respondsTo(it, "getNucleusROI")) || it.getNucleusROI() == null
}
if (!missingNucleus.isEmpty()) {
    throw new IllegalStateException(
        "${missingNucleus.size()} detection(s) have no nucleus ROI; refusing to emit an incomplete " +
        "cell/nucleus GeoJSON."
    )
}
def ids = cells.collect { it.getID()?.toString() }
if (ids.any { it == null } || ids.toSet().size() != ids.size()) {
    throw new IllegalStateException(
        "Detection UUIDs must be non-null and unique before geometry export."
    )
}

// ---- duplicate-detection guard ---------------------------------------------------------------
// Detection run with selectAnnotations() segments islet nuclei TWICE (once under Tissue, once under
// each nested Islet_N) because InstanSeg does not deduplicate across selected parents. The copies
// carry distinct UUIDs, so every uniqueness check downstream passes and the corruption stays
// invisible until it reaches the analysis. 155,669 such pairs shipped in 7 donors before this was
// caught. This is the last point where the whole detection set is in one place, so check it here.
// See scripts/groovy/0_cell_detection.groovy.
double DUP_RADIUS_UM = 1.0
double DUP_AREA_TOL  = 0.15
def cal = server.getPixelCalibration()
double pxSize = cal.hasPixelSizeMicrons() ? cal.getAveragedPixelSizeMicrons() : 1.0
double rPx = DUP_RADIUS_UM / pxSize
def buckets = [:].withDefault { [] }
for (c in cells) {
    def roi = c.getROI()
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
            if (Math.hypot(ra.getCentroidX() - rb.getCentroidX(),
                           ra.getCentroidY() - rb.getCentroidY()) > rPx) continue
            double aa = ra.getArea(), ab = rb.getArea()
            if (Math.abs(aa - ab) / Math.max(Math.max(aa, ab), 1e-9) < DUP_AREA_TOL) dupPairs++
        }
    }
}
if (dupPairs > 0) {
    throw new IllegalStateException(
        "${dupPairs} duplicate detection pair(s) found. Detection was probably run with " +
        "selectAnnotations(); re-run scripts/groovy/0_cell_detection.groovy before export."
    )
} else {
    println "[export] duplicate check: PASS"
}
// -----------------------------------------------------------------------------------------------

mkdirs(outDir)
def imageTag = (getProjectEntry()?.getImageName() ?: server.getMetadata()?.getName() ?: 'image')
        .replaceAll('[^A-Za-z0-9._-]', '_')
def outPath = buildFilePath(outDir, "cells__${imageTag}.geojson")

println "[export] ${cells.size()} cells -> ${outPath}"
long t0 = System.currentTimeMillis()

// FEATURE_COLLECTION -> {"type":"FeatureCollection","features":[...]}; feature id = UUID.
// EXCLUDE_MEASUREMENTS -> geometry + id only (small files). No GZIP, compact JSON.
exportObjectsToGeoJson(cells, outPath, "FEATURE_COLLECTION", "EXCLUDE_MEASUREMENTS")

def sizeMB = new File(outPath).length() / (1024.0 * 1024.0)
println String.format("[export] Wrote %s (%.1f MB) in %.1fs",
        outPath, sizeMB, (System.currentTimeMillis() - t0) / 1000.0)
println "[export] Done."
