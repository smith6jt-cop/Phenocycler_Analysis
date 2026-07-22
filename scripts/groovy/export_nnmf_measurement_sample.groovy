import java.nio.charset.StandardCharsets
import java.util.zip.GZIPOutputStream

def markers = [
    "CD3e", "CD20", "CD68", "CD11b", "E-cadherin",
    "CD31", "B3TUBB", "Vimentin", "EpCAM"
]
def compartments = ["Nucleus", "Cytoplasm", "Membrane", "Cell"]
def measurements = compartments.collectMany { compartment ->
    markers.collect { marker -> "${compartment}: ${marker}: Mean" }
}

def outputPath = args && args[0] ?
    args[0] :
    buildFilePath(PROJECT_BASE_DIR, "results", "nnmf_measurement_sample.tsv.gz")
def detections = getDetectionObjects()
def writer = new BufferedWriter(new OutputStreamWriter(
    new GZIPOutputStream(new FileOutputStream(outputPath)),
    StandardCharsets.UTF_8
))

try {
    writer.write((["object_id"] + measurements).join("\t"))
    writer.newLine()
    int written = 0
    for (def detection : detections) {
        def objectId = detection.getID().toString()
        if (Math.floorMod(objectId.hashCode(), 5) != 0) {
            continue
        }
        def values = measurements.collect { name ->
            def value = detection.getMeasurementList().get(name)
            Double.isFinite(value) ? Double.toString(value) : "NaN"
        }
        writer.write(([objectId] + values).join("\t"))
        writer.newLine()
        written++
    }
    println "[export] wrote ${written} sampled detections to ${outputPath}"
} finally {
    writer.close()
}
