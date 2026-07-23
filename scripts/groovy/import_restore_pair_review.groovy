import qupath.fx.dialogs.Dialogs
import qupath.lib.objects.classes.PathClass

// Template for the importer bundled with each RESTORE review package.
// Use the generated copy in expanded_review_v*/ rather than running this template.
// Import one RESTORE pair-review CSV into the currently open QuPath image.
// Existing PathClasses are backed up in object metadata before temporary review
// classes are applied, then restored exactly through the Clear dialog option.
//
// review_code:
//  -1 unavailable, 0 double-low, 1 reference control, 2 high-confidence target,
//   3 double-high, 4 other retained, 5 lower-confidence reference-negative target

// ---- config ----------------------------------------------------------------
def configuredFolder = "__RESTORE_REVIEW_CSV_DIR__"
if (configuredFolder.startsWith("__RESTORE_")) {
    println "ERROR: use the generated importer bundled with the review package"
    return
}
def folder = new File(configuredFolder)
def setTemporaryClasses = true
// ----------------------------------------------------------------------------

if (!folder.isDirectory()) {
    println "ERROR: review CSV directory not found: ${folder}"
    return
}

def currentImage = getProjectEntry() != null
    ? getProjectEntry().getImageName()
    : getCurrentImageData().getServer().getMetadata().getName()
def reviewClassPrefix = "RESTORE review - "
def originalClassMetadataKey = "RESTORE review original PathClass"
def unclassifiedSentinel = "__RESTORE_REVIEW_UNCLASSIFIED__"

def csvMetadata = { candidate ->
    def candidateReader = new BufferedReader(new FileReader(candidate))
    try {
        def candidateHeader = candidateReader.readLine()
        def candidateFirstRow = candidateReader.readLine()
        if (candidateHeader == null || candidateFirstRow == null) {
            return null
        }
        def candidateColumns = candidateHeader.split(",", -1) as List
        def imageIndex = candidateColumns.indexOf("image")
        def targetIndex = candidateColumns.indexOf("target")
        def referenceIndex = candidateColumns.indexOf("reference")
        if ([imageIndex, targetIndex, referenceIndex].any { it < 0 }) {
            return null
        }
        def candidateFields = candidateFirstRow.split(",", -1)
        if ([imageIndex, targetIndex, referenceIndex].any {
            it >= candidateFields.length
        }) {
            return null
        }
        return [
            file: candidate,
            image: candidateFields[imageIndex],
            target: candidateFields[targetIndex],
            reference: candidateFields[referenceIndex],
        ]
    } finally {
        candidateReader.close()
    }
}

def csvFiles = folder.listFiles().findAll {
    it.isFile() && it.getName().endsWith(".csv")
}
def malformedFiles = []
def allOptions = csvFiles.collect {
    def metadata = csvMetadata(it)
    if (metadata == null) {
        malformedFiles << it.getName()
    }
    return metadata
}.findAll { it != null }
if (!malformedFiles.isEmpty()) {
    println "ERROR: malformed review CSVs: ${malformedFiles.sort()}"
    return
}

def validOptions = allOptions.findAll {
    it.image == currentImage
}
if (validOptions.isEmpty()) {
    println "ERROR: this review bundle has no CSV for the open image: ${currentImage}"
    return
}

def optionByLabel = [:]
def duplicateLabels = []
validOptions.each { option ->
    def label = "${option.target} <- ${option.reference}"
    if (optionByLabel.containsKey(label)) {
        duplicateLabels << label
    } else {
        optionByLabel[label] = option
    }
}
if (!duplicateLabels.isEmpty()) {
    println "ERROR: duplicate review options for ${currentImage}: ${duplicateLabels}"
    return
}
def clearOption = "[Clear temporary RESTORE review classes]"
def choices = optionByLabel.keySet().sort() + [clearOption]
def selectedLabel = Dialogs.showChoiceDialog(
    "RESTORE pair review",
    "Open image: ${currentImage}\nChoose a target <- reference overlay:",
    choices,
    choices.first()
)
if (selectedLabel == null) {
    println "INFO: review import cancelled"
    return
}

if (selectedLabel == clearOption) {
    int restoredCount = 0
    int legacyClearedCount = 0
    getDetectionObjects().each { object ->
        def metadata = object.getMetadata()
        if (metadata.containsKey(originalClassMetadataKey)) {
            def originalClass = metadata.remove(originalClassMetadataKey)
            object.setPathClass(
                originalClass == unclassifiedSentinel
                    ? null
                    : PathClass.fromString(originalClass)
            )
            restoredCount++
        } else if (object.getPathClass() != null &&
                   object.getPathClass().toString().startsWith(reviewClassPrefix)) {
            // Importers predating the metadata backup only replaced unclassified cells.
            object.setPathClass(null)
            legacyClearedCount++
        }
    }
    fireHierarchyUpdate()
    println "Restored original PathClasses for ${restoredCount} detections."
    if (legacyClearedCount > 0) {
        println "Cleared ${legacyClearedCount} legacy RESTORE review classes to Unclassified."
    }
    return
}

def file = optionByLabel[selectedLabel].file
println "INFO: selected ${selectedLabel} (${file.getName()})"

def classColors = [
    "RESTORE review - Unavailable"       : [34, 34, 34],
    "RESTORE review - Double-low"        : [189, 189, 189],
    "RESTORE review - Reference control" : [0, 114, 178],
    "RESTORE review - Called target+"          : [213, 94, 0],
    "RESTORE review - Retained below Step 2"   : [185, 175, 203],
    "RESTORE review - Double-high"       : [204, 121, 167],
    "RESTORE review - Other retained"    : [127, 127, 127],
]
def codeClasses = [
    "-1": "RESTORE review - Unavailable",
    "0" : "RESTORE review - Double-low",
    "1" : "RESTORE review - Reference control",
    "2" : "RESTORE review - Called target+",
    "3" : "RESTORE review - Double-high",
    "4" : "RESTORE review - Other retained",
    "5" : "RESTORE review - Retained below Step 2",
]
def hasTemporaryReviewClass = { object ->
    def pathClass = object.getPathClass()
    pathClass != null && pathClass.toString().startsWith(reviewClassPrefix)
}

def reader = new BufferedReader(new FileReader(file))
def headerLine = reader.readLine()
if (headerLine == null) {
    reader.close()
    println "ERROR: review CSV is empty: ${file}"
    return
}
def header = headerLine.split(",", -1)
def column = [:]
header.eachWithIndex { name, index -> column[name] = index }
def required = ["object_id", "image", "review_code", "target", "reference"]
if (!required.every { column.containsKey(it) }) {
    reader.close()
    println "ERROR: missing required columns; found ${header as List}"
    return
}

def rows = []
String line
while ((line = reader.readLine()) != null) {
    def fields = line.split(",", -1)
    if (fields.length != header.length) {
        reader.close()
        println "ERROR: malformed CSV row with ${fields.length} fields; expected ${header.length}"
        return
    }
    rows << fields
}
reader.close()
if (rows.isEmpty()) {
    println "ERROR: review CSV has no data rows: ${file}"
    return
}

def expectedImages = rows.collect { it[column["image"]] }.toSet()
if (expectedImages.size() != 1) {
    println "ERROR: review CSV contains ${expectedImages.size()} image names; expected exactly one"
    return
}
def expectedImage = expectedImages.first()
if (expectedImage != currentImage) {
    println "ERROR: CSV/image mismatch; nothing was imported."
    println "  CSV expects: ${expectedImage}"
    println "  Open image:  ${currentImage}"
    if (!validOptions.isEmpty()) {
        println "Valid review CSVs for the open image:"
        validOptions.each {
            println "  ${it.target} <- ${it.reference} (${it.file.getName()})"
        }
    }
    return
}

def targetValues = rows.collect { it[column["target"]] }.toSet()
def referenceValues = rows.collect { it[column["reference"]] }.toSet()
if (targetValues.size() != 1 || referenceValues.size() != 1) {
    println "ERROR: review CSV must contain exactly one target/reference pair"
    return
}
def target = targetValues.first()
def reference = referenceValues.first()
def measurementName = "RESTORE review: ${target} <- ${reference}"

def detections = [:]
getDetectionObjects().each { detections[it.getID().toString()] = it }
if (detections.isEmpty()) {
    println "ERROR: the current image has no detections"
    return
}

def resolved = []
def missingIds = []
for (def fields : rows) {
    def objectId = fields[column["object_id"]]
    def object = detections[objectId]
    if (object == null) {
        if (missingIds.size() < 5) {
            missingIds << objectId
        }
        continue
    }
    def codeText = fields[column["review_code"]]
    if (!codeClasses.containsKey(codeText)) {
        println "ERROR: unsupported review_code '${codeText}'; nothing was imported"
        return
    }
    resolved << [objectId: objectId, object: object, code: codeText]
}

if (resolved.size() != rows.size()) {
    println "ERROR: UUID mismatch within the correct image; nothing was imported."
    println "  Matched ${resolved.size()} of ${rows.size()} review rows."
    println "  Example missing IDs: ${missingIds}"
    println "The image detections may have been regenerated after the review data were exported."
    return
}

if (setTemporaryClasses) {
    classColors.each { name, color ->
        getPathClass(name).setColor(color[0], color[1], color[2])
    }
}

int backedUpCount = 0
int legacyReviewCount = 0
if (setTemporaryClasses) {
    def objectsToClassify = resolved.collect { it.object }
    objectsToClassify.each { object ->
        def metadata = object.getMetadata()
        if (!metadata.containsKey(originalClassMetadataKey)) {
            def pathClass = object.getPathClass()
            if (hasTemporaryReviewClass(object)) {
                // Earlier importer versions only overlaid cells that were unclassified.
                metadata[originalClassMetadataKey] = unclassifiedSentinel
                legacyReviewCount++
            } else {
                metadata[originalClassMetadataKey] = pathClass == null
                    ? unclassifiedSentinel
                    : pathClass.toString()
            }
            backedUpCount++
        }
    }
}

resolved.each { match ->
    match.object.getMeasurementList().put(measurementName, match.code as double)
    if (setTemporaryClasses) {
        match.object.setPathClass(getPathClass(codeClasses[match.code]))
    }
}
fireHierarchyUpdate()

println "Imported '${measurementName}' for ${resolved.size()} detections (0 unmatched)."
if (setTemporaryClasses) {
    println "Applied temporary RESTORE review PathClasses."
    println "Backed up ${backedUpCount} original PathClasses for exact restoration."
    if (legacyReviewCount > 0) {
        println "Migrated ${legacyReviewCount} review classes from the legacy importer."
    }
    println "Cells outside the representative crop retained their existing PathClasses."
    println "Choose '${clearOption}' when finished to restore original PathClasses."
} else {
    println "PathClasses were intentionally preserved. Color detections by '${measurementName}'."
}
