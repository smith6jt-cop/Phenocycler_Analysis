// export_cells_geojson.groovy  (QuPath 0.5+)
// Export full-resolution per-cell boundaries as GeoJSON for REDSEA (phenocycler.redsea).
//
// REDSEA rasterizes these polygons into an int32 instance mask and joins the corrected
// per-cell means back to QuPath's detection UUID, so the GeoJSON MUST:
//   * be in FULL-RESOLUTION PIXEL coordinates (not microns, no downsample), and
//   * carry each feature's id == the detection's PathObject.getID() (a UUID).
// QuPath's GeoJSON exporter writes `PathObject.getID()` as the feature `id` and pixel
// coordinates by default, which is exactly what redsea_full.resolve_paths expects.
//
// Output name must match redsea.resolve_paths: cells__<tag>.geojson where <tag> is the
// QuPath image name with [^A-Za-z0-9._-] -> '_'  (e.g. "6539_Scan1.er.qptiff - resolution #1").
//
// Run headlessly per image (GUI closed):
//   QuPath script --project <project.qpproj> --image "<image name>" export_cells_geojson.groovy
// See scripts/export_new_geojsons.sh for the batch driver.

import qupath.lib.io.GsonTools
import static qupath.lib.scripting.QP.*

// ---- config ----------------------------------------------------------------
// <data_dir>/redsea_scratch/geojson (edit to match config.ini [paths] data_dir)
def outDir = new File(System.getProperty("user.home"),
                      "Phenocycler_Analysis/data/redsea_scratch/geojson")
// ----------------------------------------------------------------------------

outDir.mkdirs()
def imageName = getProjectEntry() != null ? getProjectEntry().getImageName()
                                          : getCurrentImageData().getServer().getMetadata().getName()
def tag = imageName.replaceAll("[^A-Za-z0-9._-]", "_")
def outFile = new File(outDir, "cells__${tag}.geojson")

def detections = getDetectionObjects()
if (detections.isEmpty()) { println "ERROR: no detections in ${imageName} — run cell detection first"; return }

// includeMeasurements=false keeps the file small; FEATURE_COLLECTION + pretty for readability.
def gson = GsonTools.getInstance(true)
outFile.withWriter('UTF-8') { w ->
    exportObjectsToGeoJson(detections, outFile.getAbsolutePath(),
        "FEATURE_COLLECTION", "EXCLUDE_MEASUREMENTS")
}
println "Wrote ${detections.size()} cell polygons -> ${outFile}"
