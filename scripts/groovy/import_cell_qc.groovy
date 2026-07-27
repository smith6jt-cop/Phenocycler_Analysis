import qupath.fx.dialogs.Dialogs
import qupath.lib.objects.classes.PathClass

// Highlight the detections that single-cell QC removed, in the currently open QuPath image.
//
// Matching is on the QuPath detection UUID (PathObject.getID()), the same key REDSEA, RESTORE and
// the lineage importer use — so there is no centroid rounding and the match is exact.
//
// Existing PathClasses are backed up into object metadata BEFORE any temporary class is applied,
// and the Clear option restores every one of them exactly, including Unclassified. Running this
// script never destroys existing classification work.
//
// qc_code:
//   0 retained (only present if the CSV was exported with --include-retained)
//   1 no nucleus          - QuPath drew a cell but found no nucleus in it (~1/3 normal DNA stain)
//   2 area below min      - smaller than the configured floor; fragments, dim in DNA
//   3 area above max      - > 4x this donor's median; merged stromal / vessel-wall blobs
//   4 low solidity        - deeply concave or multi-lobed; usually two cells traced as one
//   5 raster degenerate   - outline claims one area but fills to a small fraction of it
//   6 raster empty        - outline fills to zero pixels
//   7 no cytoplasm        - RETAINED and still called; nucleus == cell, so it may not set a
//                           threshold. Shown so the population is visible, not because it is bad.

// ---- config ----------------------------------------------------------------
def configuredFolder = "__CELL_QC_CSV_DIR__"
if (configuredFolder.startsWith("__CELL_QC_")) {
    println "ERROR: run the generated copy written by `python -m phenocycler.cell_qc_export`,"
    println "       not this template (it has no CSV directory baked in)."
    return
}
def folder = new File(configuredFolder)
def setTemporaryClasses = true
// ----------------------------------------------------------------------------

if (!folder.isDirectory()) {
    println "ERROR: cell-QC CSV directory not found: ${folder}"
    return
}

def qcClassPrefix = "Cell QC - "
def originalClassMetadataKey = "Cell QC original PathClass"
def unclassifiedSentinel = "__CELL_QC_UNCLASSIFIED__"

def codeClasses = [
    0: qcClassPrefix + "Retained",
    1: qcClassPrefix + "No nucleus",
    2: qcClassPrefix + "Area below min",
    3: qcClassPrefix + "Area above max",
    4: qcClassPrefix + "Low solidity",
    5: qcClassPrefix + "Raster degenerate",
    6: qcClassPrefix + "Raster empty",
    7: qcClassPrefix + "No cytoplasm (kept)",
]
// Paul Tol bright — colourblind-safe, and distinct from the RESTORE review palette.
// RGB triples with setColor(r, g, b), matching import_broad_lineage.groovy and
// import_restore_pair_review.groovy: PathClass is immutable in QuPath 0.4+, so the colour is set on
// the interned instance returned by getPathClass(), not on a per-object copy.
def codeColors = [
    0: [190, 190, 190],   // grey    - retained background
    1: [204,  51, 119],   // magenta - no nucleus (the largest group)
    2: [238, 119,  51],   // orange  - too small
    3: [  0, 119, 187],   // blue    - too large
    4: [238, 102, 119],   // red     - bad shape
    5: [  0, 153, 136],   // teal    - broken outline
    6: [ 51,  34, 136],   // indigo  - empty outline
    7: [221, 204, 119],   // sand    - kept, estimation-excluded
]
codeColors.each { code, c -> getPathClass(codeClasses[code]).setColor(c[0], c[1], c[2]) }

def currentImage = getProjectEntry() != null
    ? getProjectEntry().getImageName()
    : getCurrentImageData().getServer().getMetadata().getName()

def detections = getDetectionObjects()
if (detections.isEmpty()) {
    println "ERROR: no detections in the open image — is this the right image?"
    return
}

// ---- restore-original-classes path -----------------------------------------
def clearOption = "[Clear temporary Cell QC classes]"

def clearTemporaryClasses = {
    int restored = 0
    int legacy = 0
    for (object in detections) {
        def metadata = object.getMetadata()
        if (metadata.containsKey(originalClassMetadataKey)) {
            def original = metadata.remove(originalClassMetadataKey)
            object.setPathClass(
                original == unclassifiedSentinel ? null : getPathClass(original)
            )
            restored++
        } else if (object.getPathClass() != null
                   && object.getPathClass().toString().startsWith(qcClassPrefix)) {
            object.setPathClass(null)
            legacy++
        }
    }
    fireHierarchyUpdate()
    println "Restored original PathClasses for ${restored} detections."
    if (legacy > 0) {
        println "Cleared ${legacy} legacy Cell QC classes to Unclassified."
    }
}

// ---- find the CSV for the open image ---------------------------------------
def csvFiles = folder.listFiles({ f ->
    f.getName().startsWith("cell_qc_") && f.getName().endsWith(".csv")
} as FileFilter)
if (csvFiles == null || csvFiles.length == 0) {
    println "ERROR: no cell_qc_*.csv found in ${folder}"
    return
}

def matching = []
for (candidate in csvFiles) {
    def reader = new BufferedReader(new FileReader(candidate))
    try {
        def header = reader.readLine()
        def firstRow = reader.readLine()
        if (header == null || firstRow == null) continue
        def columns = header.split(",", -1) as List
        def imageIndex = columns.indexOf("image")
        if (imageIndex < 0) continue
        def fields = firstRow.split(",", -1)
        if (imageIndex >= fields.length) continue
        if (fields[imageIndex] == currentImage) {
            matching << candidate
        }
    } finally {
        reader.close()
    }
}

if (matching.isEmpty()) {
    println "ERROR: no cell-QC CSV matches the open image."
    println "       open image : ${currentImage}"
    println "       available  : " + csvFiles.collect { it.getName() }.join(", ")
    println "       (the CSV's 'image' column must equal the QuPath image name exactly)"
    return
}

def choices = [clearOption] + matching.collect { it.getName() }.sort()
def chosen = Dialogs.showChoiceDialog(
    "Cell QC overlay",
    "Image: ${currentImage}\n\nApply a QC overlay, or clear a previous one.",
    choices as String[],
    choices[0]
)
if (chosen == null) {
    println "Cancelled."
    return
}
if (chosen == clearOption) {
    clearTemporaryClasses()
    return
}

def csvFile = matching.find { it.getName() == chosen }

// ---- read the CSV ----------------------------------------------------------
def codeById = [:]
def reasonCounts = [:]
def reader = new BufferedReader(new FileReader(csvFile))
int rows = 0
try {
    def columns = reader.readLine().split(",", -1) as List
    def idIndex = columns.indexOf("object_id")
    def codeIndex = columns.indexOf("qc_code")
    def reasonIndex = columns.indexOf("qc_reason")
    if (idIndex < 0 || codeIndex < 0) {
        println "ERROR: ${csvFile.getName()} is missing object_id or qc_code"
        return
    }
    def line
    while ((line = reader.readLine()) != null) {
        if (line.trim().isEmpty()) continue
        def f = line.split(",", -1)
        if (idIndex >= f.length || codeIndex >= f.length) continue
        codeById[f[idIndex]] = f[codeIndex] as int
        if (reasonIndex >= 0 && reasonIndex < f.length) {
            reasonCounts[f[reasonIndex]] = (reasonCounts[f[reasonIndex]] ?: 0) + 1
        }
        rows++
    }
} finally {
    reader.close()
}
println "Read ${rows} QC rows from ${csvFile.getName()}"

// ---- match on detection UUID ------------------------------------------------
def matches = []
int unmatchedInImage = 0
for (object in detections) {
    def code = codeById.get(object.getID().toString())
    if (code == null) {
        unmatchedInImage++
    } else {
        matches << [object: object, code: code]
    }
}

println "Matched ${matches.size()} of ${detections.size()} detections in this image."
if (matches.isEmpty()) {
    println "ERROR: no detection UUID matched. This CSV is for a different image or a different"
    println "       QuPath project — refusing to modify anything."
    return
}
int unmatchedInCsv = codeById.size() - matches.size()
if (unmatchedInCsv > 0) {
    println "NOTE: ${unmatchedInCsv} CSV rows have no detection in this image."
    println "      Expected if the project was re-segmented after the QC export."
}

println "\nQC breakdown for this image:"
reasonCounts.sort { -it.value }.each { reason, count ->
    println String.format("  %-20s %,10d", reason, count)
}
if (unmatchedInImage > 0) {
    println String.format("  %-20s %,10d", "(not in CSV)", unmatchedInImage)
}

if (!setTemporaryClasses) {
    println "\nsetTemporaryClasses = false — reported only, nothing modified."
    return
}

// ---- back up originals, then apply -----------------------------------------
for (object in detections) {
    def metadata = object.getMetadata()
    if (!metadata.containsKey(originalClassMetadataKey)) {
        def pathClass = object.getPathClass()
        if (pathClass == null) {
            metadata[originalClassMetadataKey] = unclassifiedSentinel
        } else {
            metadata[originalClassMetadataKey] = pathClass.toString()
        }
    }
}

int applied = 0
for (match in matches) {
    def name = codeClasses[match.code]
    if (name == null) continue
    match.object.setPathClass(getPathClass(name))
    applied++
}
fireHierarchyUpdate()

println "\nApplied ${applied} temporary QC classes."
println "Detections NOT in the CSV keep their original class — with an excluded-only export that"
println "means everything unhighlighted passed QC."
println "Re-run this script and choose '${clearOption}' to restore every original PathClass."
