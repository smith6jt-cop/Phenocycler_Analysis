/*
 * import_cell_assignments.groovy  (QuPath 0.5+)
 *
 * Import the uncertainty-preserving CSV written by phenocycler.qupath_export.
 * Matching uses the QuPath detection UUID (`object_id`) and the exact image
 * name. Before changing any PathClass, the script validates every exported
 * canonical UUID against the open image. QuPath detections below the ingest
 * structural floor are intentionally outside that canonical universe and are
 * reported but left unchanged.
 *
 * Required columns:
 *   object_id, image, broad_type, specific_type, assignment_status,
 *   confidence, best_broad_type, best_specific_type,
 *   broad_assignment_status, specific_assignment_status,
 *   ranked_broad_probabilities, ranked_specific_probabilities
 *
 * At the requested import level, `anchor`, `inferred`, and `other`
 * assignments receive a biological class. `ambiguous` and `unavailable`
 * NEVER receive the corresponding best candidate; they receive explicit
 * uncertainty classes. Broad and specific decisions/candidates remain
 * separate metadata, so subtype uncertainty cannot rewrite an authoritative
 * broad assignment when importing at the broad level.
 *
 * Run with the immutable run's export directory as the first script argument:
 *   QuPath script --project project.qpproj --image "exact image name" \
 *     scripts/groovy/import_cell_assignments.groovy \
 *     --args "/path/to/data/runs/<run_id>/qupath_export,specific"
 *
 * The optional second argument is `specific` (default) or `broad`.
 */

import groovy.json.JsonSlurper
import java.security.MessageDigest

// ---- configuration ---------------------------------------------------------
String assignmentFolder = (
    args != null && args.size() > 0 ? args[0]?.trim() : ""
)
String assignmentLevel = (
    args != null && args.size() > 1 ? args[1]?.trim()?.toLowerCase() : "specific"
)
boolean storeRankedProbabilities = false
// ----------------------------------------------------------------------------

if (assignmentFolder == null || assignmentFolder.isEmpty()) {
    println "ERROR: pass the run's qupath_export directory as the first script argument."
    return
}
if (!(assignmentLevel in ["broad", "specific"])) {
    println "ERROR: assignment level must be 'broad' or 'specific', got '${assignmentLevel}'."
    return
}
def folder = new File(assignmentFolder)
if (!folder.isDirectory()) {
    println "ERROR: assignment directory not found: ${folder}"
    return
}

def requiredColumns = [
    "object_id", "image", "broad_type", "specific_type",
    "assignment_status", "confidence",
    "best_broad_type", "best_specific_type",
    "broad_assignment_status", "specific_assignment_status",
    "ranked_broad_probabilities", "ranked_specific_probabilities",
    "ranked_probabilities"
]
def acceptedStatuses = ["anchor", "inferred", "other"] as Set
def uncertaintyStatuses = ["ambiguous", "unavailable"] as Set

// RFC-4180 field parsing, including the quoted JSON in ranked_probabilities.
// Pandas does not emit embedded newlines for this export, so parsing is safely
// line-oriented.
def parseCsvLine = { String line ->
    def fields = []
    def current = new StringBuilder()
    boolean quoted = false
    for (int i = 0; i < line.length(); i++) {
        char c = line.charAt(i)
        if (((int) c) == 34) {
            if (quoted && i + 1 < line.length() && ((int) line.charAt(i + 1)) == 34) {
                current.append((char) 34)
                i++
            } else {
                quoted = !quoted
            }
        } else if (((int) c) == 44 && !quoted) {
            fields << current.toString()
            current.setLength(0)
        } else {
            current.append(c)
        }
    }
    if (quoted) {
        throw new IllegalArgumentException("unterminated quoted CSV field")
    }
    fields << current.toString()
    return fields
}

def currentImage = getProjectEntry() != null
    ? getProjectEntry().getImageName()
    : getCurrentImageData().getServer().getMetadata().getName()
def detections = getDetectionObjects()
if (detections.isEmpty()) {
    println "ERROR: no detections in the open image."
    return
}
def detectionsById = [:]
for (object in detections) {
    def id = object.getID()?.toString()
    if (id == null || detectionsById.containsKey(id)) {
        println "ERROR: detection UUIDs must be non-null and unique."
        return
    }
    detectionsById[id] = object
}

def csvFiles = folder.listFiles({ file ->
    file.isFile() &&
        file.getName().startsWith("cell_assignments_") &&
        file.getName().endsWith(".csv")
} as FileFilter)
if (csvFiles == null || csvFiles.length == 0) {
    println "ERROR: no cell_assignments_*.csv files in ${folder}"
    return
}

// Identify the file from its authoritative image column, not its filename.
def matchingFiles = []
for (file in csvFiles) {
    def reader = new BufferedReader(new FileReader(file))
    try {
        def headerLine = reader.readLine()
        def firstLine = reader.readLine()
        if (headerLine == null || firstLine == null) continue
        def header = parseCsvLine(headerLine)
        int imageIndex = header.indexOf("image")
        if (imageIndex < 0) continue
        def first = parseCsvLine(firstLine)
        if (imageIndex < first.size() && first[imageIndex] == currentImage) {
            matchingFiles << file
        }
    } finally {
        reader.close()
    }
}
if (matchingFiles.size() != 1) {
    println "ERROR: expected exactly one assignment CSV for '${currentImage}', found ${matchingFiles.size()}."
    println "       candidates: " + matchingFiles.collect { it.getName() }.join(", ")
    return
}
def csvFile = matchingFiles[0]

// Verify the exact immutable export before interpreting any row. This catches
// a truncated or manually edited CSV even when every remaining UUID is valid.
def exportManifestFile = new File(
    folder.getParentFile(), "manifests/qupath_export.json"
)
if (!exportManifestFile.isFile()) {
    println "ERROR: export manifest not found: ${exportManifestFile}"
    return
}
def exportManifest
try {
    exportManifest = new JsonSlurper().parse(exportManifestFile)
} catch (Exception error) {
    println "ERROR: cannot parse export manifest: ${error.getMessage()}"
    return
}
def fingerprintMatches = exportManifest.files.findAll { item ->
    new File(item.path.toString()).getCanonicalFile() ==
        csvFile.getCanonicalFile()
}
if (fingerprintMatches.size() != 1) {
    println "ERROR: export manifest does not identify ${csvFile.getCanonicalPath()} exactly once."
    return
}
def fingerprint = fingerprintMatches[0]
if (fingerprint.hash_mode != "full") {
    println "ERROR: assignment CSV fingerprint must use a full SHA-256 hash."
    return
}
if (csvFile.length() != (fingerprint.size_bytes as long)) {
    println "ERROR: assignment CSV size differs from its immutable manifest."
    return
}
def digest = MessageDigest.getInstance("SHA-256")
csvFile.withInputStream { stream ->
    byte[] buffer = new byte[1024 * 1024]
    int count
    while ((count = stream.read(buffer)) > 0) {
        digest.update(buffer, 0, count)
    }
}
def actualSha256 = digest.digest().collect {
    String.format("%02x", it & 0xff)
}.join()
if (actualSha256 != fingerprint.sha256) {
    println "ERROR: assignment CSV SHA-256 differs from its immutable manifest."
    return
}

// Reusable streaming reader. The first pass validates everything; the second
// applies rows only after the contract is proven.
def visitRows = { Closure visitor ->
    def reader = new BufferedReader(new FileReader(csvFile))
    try {
        def headerLine = reader.readLine()
        if (headerLine == null) {
            throw new IllegalArgumentException("empty assignment CSV")
        }
        def header = parseCsvLine(headerLine)
        if (header.toSet().size() != header.size()) {
            throw new IllegalArgumentException("duplicate CSV column names")
        }
        def missing = requiredColumns.findAll { !header.contains(it) }
        if (!missing.isEmpty()) {
            throw new IllegalArgumentException("missing columns: ${missing}")
        }
        def column = [:]
        header.eachWithIndex { name, index -> column[name] = index }
        String line
        int lineNumber = 1
        while ((line = reader.readLine()) != null) {
            lineNumber++
            if (line.trim().isEmpty()) continue
            def fields = parseCsvLine(line)
            if (fields.size() != header.size()) {
                throw new IllegalArgumentException(
                    "line ${lineNumber}: ${fields.size()} fields, expected ${header.size()}"
                )
            }
            visitor(fields, column, lineNumber)
        }
    } finally {
        reader.close()
    }
}

def seen = new HashSet<String>()
int rowCount = 0
try {
    visitRows { fields, column, lineNumber ->
        def id = fields[column["object_id"]]
        def image = fields[column["image"]]
        def statusColumn = (
            assignmentLevel == "broad"
                ? "broad_assignment_status"
                : "specific_assignment_status"
        )
        def status = fields[column[statusColumn]]
        if (image != currentImage) {
            throw new IllegalArgumentException(
                "line ${lineNumber}: image '${image}' != open image '${currentImage}'"
            )
        }
        if (!detectionsById.containsKey(id)) {
            throw new IllegalArgumentException(
                "line ${lineNumber}: object_id ${id} is absent from the open image"
            )
        }
        if (!seen.add(id)) {
            throw new IllegalArgumentException(
                "line ${lineNumber}: duplicate object_id ${id}"
            )
        }
        if (!(status in acceptedStatuses) && !(status in uncertaintyStatuses)) {
            throw new IllegalArgumentException(
                "line ${lineNumber}: unsupported assignment_status '${status}'"
            )
        }
        def broad = fields[column["broad_type"]]
        def specific = fields[column["specific_type"]]
        if (status in acceptedStatuses && (broad == null || broad.isEmpty())) {
            throw new IllegalArgumentException(
                "line ${lineNumber}: accepted assignment has no broad_type"
            )
        }
        if (
            status in acceptedStatuses &&
            assignmentLevel == "specific" &&
            (specific == null || specific.isEmpty())
        ) {
            throw new IllegalArgumentException(
                "line ${lineNumber}: accepted assignment has no specific_type"
            )
        }
        def confidenceText = fields[column["confidence"]]
        if (confidenceText != null && !confidenceText.isEmpty()) {
            double confidence = Double.parseDouble(confidenceText)
            if (!Double.isFinite(confidence) || confidence < 0.0 || confidence > 1.0) {
                throw new IllegalArgumentException(
                    "line ${lineNumber}: confidence must be in [0,1]"
                )
            }
        }
        rowCount++
    }
} catch (Exception error) {
    println "ERROR: assignment validation failed: ${error.getMessage()}"
    return
}
if (seen.isEmpty()) {
    println "ERROR: assignment CSV contains no canonical detections."
    return
}
def outsideCanonicalUniverse = detectionsById.keySet().findAll {
    !seen.contains(it)
}
println "Validated ${rowCount} exact canonical UUID assignments for '${currentImage}'."
if (!outsideCanonicalUniverse.isEmpty()) {
    println "Leaving ${outsideCanonicalUniverse.size()} non-canonical QuPath detections unchanged " +
        "(for example, objects below the ingest structural floor)."
}

def broadColors = [
    "Immune"      : [34, 136, 51],
    "Endothelial" : [102, 204, 238],
    "Epithelial"  : [68, 119, 170],
    "Neural"      : [181, 131, 141],
    "Mesenchymal" : [170, 51, 119],
    "Other"       : [187, 187, 187],
]
def uncertaintyColors = [
    "Phenocycler - Ambiguous"   : [238, 119, 51],
    "Phenocycler - Unavailable" : [102, 102, 102],
]
def previousClassKey = "Phenocycler previous PathClass"
def metadataPrefix = "Phenocycler "
int applied = 0

visitRows { fields, column, lineNumber ->
    def id = fields[column["object_id"]]
    def broad = fields[column["broad_type"]]
    def specific = fields[column["specific_type"]]
    def combinedStatus = fields[column["assignment_status"]]
    def broadStatus = fields[column["broad_assignment_status"]]
    def specificStatus = fields[column["specific_assignment_status"]]
    def status = assignmentLevel == "broad" ? broadStatus : specificStatus
    def bestBroad = fields[column["best_broad_type"]]
    def bestSpecific = fields[column["best_specific_type"]]
    def broadRanking = fields[column["ranked_broad_probabilities"]]
    def specificRanking = fields[column["ranked_specific_probabilities"]]
    def confidenceText = fields[column["confidence"]]
    double confidence = (
        confidenceText == null || confidenceText.isEmpty()
            ? Double.NaN
            : Double.parseDouble(confidenceText)
    )

    String className
    def color
    if (status == "ambiguous") {
        className = "Phenocycler - Ambiguous"
        color = uncertaintyColors[className]
    } else if (status == "unavailable") {
        className = "Phenocycler - Unavailable"
        color = uncertaintyColors[className]
    } else {
        if (broad == null || broad.isEmpty()) {
            throw new IllegalArgumentException(
                "line ${lineNumber}: accepted assignment has no broad_type"
            )
        }
        className = broad
        if (
            assignmentLevel == "specific" &&
            specific != null &&
            !specific.isEmpty() &&
            specific != broad
        ) {
            className = "${broad} - ${specific}"
        }
        color = broadColors[broad] ?: [120, 120, 120]
    }

    def pathClass = getPathClass(className)
    pathClass.setColor(color[0], color[1], color[2])
    def object = detectionsById[id]
    def metadata = object.getMetadata()
    if (!metadata.containsKey(previousClassKey)) {
        metadata[previousClassKey] = (
            object.getPathClass() == null
                ? "__UNCLASSIFIED__"
                : object.getPathClass().toString()
        )
    }
    object.setPathClass(pathClass)
    object.getMeasurementList().put("Phenocycler confidence", confidence)
    metadata[metadataPrefix + "assignment level"] = assignmentLevel
    metadata[metadataPrefix + "assignment status"] = combinedStatus
    metadata[metadataPrefix + "broad assignment status"] = broadStatus
    metadata[metadataPrefix + "specific assignment status"] = specificStatus
    metadata[metadataPrefix + "broad type"] = broad
    metadata[metadataPrefix + "specific type"] = specific
    metadata[metadataPrefix + "best broad candidate"] = bestBroad
    metadata[metadataPrefix + "best specific candidate"] = bestSpecific
    if (storeRankedProbabilities) {
        metadata[metadataPrefix + "ranked broad probabilities"] = broadRanking
        metadata[metadataPrefix + "ranked specific probabilities"] = specificRanking
    }
    applied++
}

fireHierarchyUpdate()
println "Imported ${applied} assignments from ${csvFile.getName()} at '${assignmentLevel}' level."
println "Ambiguous/unavailable detections retain their candidates as metadata but use uncertainty classes."
