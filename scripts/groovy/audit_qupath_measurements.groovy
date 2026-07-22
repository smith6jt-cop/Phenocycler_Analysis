def detections = getDetectionObjects()
println "[audit] detections=${detections.size()}"
if (detections.isEmpty()) {
    return
}

def wantedMarkers = [
    "CD3e", "CD20", "CD68", "CD11b", "E-cadherin",
    "CD31", "B3TUBB", "Vimentin", "EpCAM"
] as Set
def wantedCompartments = ["Nucleus", "Cytoplasm", "Membrane", "Cell"] as Set
def names = detections[0].getMeasurementList().getNames() as List
def selected = names.findAll { name ->
    def parts = name.split(": ", 3)
    parts.size() == 3 &&
        wantedCompartments.contains(parts[0]) &&
        wantedMarkers.contains(parts[1])
}

println "[audit] total_measurements=${names.size()}"
println "[audit] selected_measurements=${selected.size()}"
selected.sort().each { println "[measurement] ${it}" }

for (int i = 0; i < Math.min(3, detections.size()); i++) {
    def detection = detections[i]
    println "[cell] index=${i} id=${detection.getID()}"
    selected.findAll { it.endsWith(": Mean") }.sort().each { name ->
        println "[value] ${name}=${detection.getMeasurementList().get(name)}"
    }
}
