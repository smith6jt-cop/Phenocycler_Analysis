/*
 * export_cells_geojson.groovy
 *
 * Per-CELL geometry export for the scalable pixel-level REDSEA pipeline
 * (scripts/senior/redsea_full.py). Adapted from
 * export_isletsAndexpansions_geojson.groovy: same exporter, same FeatureCollection
 * format, but exports the cell DETECTIONS (not Islet annotations) and EXCLUDES
 * measurements (we already have those in data/cells/*.parquet — we only need the
 * polygon geometry + the object UUID for the join).
 *
 * Each exported feature carries:
 *   - "id"        : the QuPath detection UUID == object_id in data/cells/*.parquet
 *   - "geometry"  : the cell polygon, in FULL-RES pixel coordinates (downsample 1)
 *
 * Output (one file per image):
 *   <OUT_DIR>/cells__<imageTag>.geojson
 *
 * Run headless, fresh JVM (QuPath GUI must be closed so the heap is free), e.g.:
 *   /home/smith6jt/QuPath/bin/QuPath script \
 *       --project /home/smith6jt/IO60panc2nd/project.qpproj \
 *       --image "6539_Scan1.er.qptiff - resolution #1" \
 *       /home/smith6jt/IO60panc2nd/scripts/export_cells_geojson.groovy
 *
 * Optional: pass an output dir as the first script arg (--args "/some/dir").
 */

// ==================== PARAMETERS ====================
String DEFAULT_OUT_DIR =
    '/home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/data/redsea_scratch/geojson'
// ====================================================

String outDir = (args != null && args.size() > 0 && args[0]?.trim()) ? args[0].trim()
                                                                      : DEFAULT_OUT_DIR

def server = getCurrentServer()
if (server == null) { println "ERROR: no image open."; return }

// Cells = detection objects with a ROI. (Annotations like Islet/Tissue are skipped.)
def cells = getDetectionObjects().findAll { it.getROI() != null }
if (cells.isEmpty()) {
    println "ERROR: no detection objects (cells) found in this image."
    return
}

// ---- duplicate-detection guard ---------------------------------------------------------------
// Detection run with selectAnnotations() segments islet nuclei TWICE (once under Tissue, once under
// each nested Islet_N) because InstanSeg does not deduplicate across selected parents. The copies
// carry distinct UUIDs, so every uniqueness check downstream passes and the corruption stays
// invisible until it reaches the analysis. 155,669 such pairs shipped in 7 donors before this was
// caught. This is the last point where the whole detection set is in one place, so check it here.
// See scripts/0_cell_detection.groovy and phenocycler.cell_qc.find_duplicate_detections.
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
    println "[export] WARNING: ${dupPairs} duplicate detection pair(s) in this image."
    println "[export]          Detection was probably run with selectAnnotations() instead of"
    println "[export]          selectObjectsByClassification(\"Tissue\"). Re-run detection with"
    println "[export]          scripts/0_cell_detection.groovy, or the duplicates get exported."
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
