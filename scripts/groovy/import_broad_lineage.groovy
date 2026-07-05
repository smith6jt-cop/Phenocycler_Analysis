// import_broad_lineage.groovy  (QuPath 0.5+)
// Set each detection's PathClass to its broad lineage, matched by UUID.
//
// QuPath's PathObject.getID() == our object_id, so matching is EXACT (no centroid
// rounding) and the class is set in ONE step.
//
// CSV format (phenocycler.qupath_export):  object_id, broad_lineage, image
//
// EDIT `folder` below to point at <data_dir>/phenotype/qupath_class produced by the
// pipeline (python -m phenocycler.qupath_export).

// ---- config ----------------------------------------------------------------
def donorsToImport = []          // [] / null = every pheno_class_*.csv in the folder
def folder = new File(System.getProperty("user.home"),
                      "Phenocycler_Analysis/data/phenotype/qupath_class")
// ----------------------------------------------------------------------------

// distinct colours for the six broad lineages (match the analysis figures)
def COLORS = [
    "Epithelial" : [68, 119, 170],   "Fibroblast" : [238, 102, 119],
    "Immune"     : [34, 136, 51],    "Endocrine"  : [204, 187, 68],
    "Endothelial": [102, 204, 238],  "Muscle"     : [170, 51, 119],
]
COLORS.each { name, c -> getPathClass(name).setColor(c[0], c[1], c[2]) }

if (!folder.isDirectory()) { println "ERROR: folder not found: ${folder}"; return }
def csvFiles = folder.listFiles().findAll {
    it.isFile() && it.getName().startsWith("pheno_class_") && it.getName().endsWith(".csv") &&
    (donorsToImport == null || donorsToImport.isEmpty() ||
     donorsToImport.any { d -> it.getName() == "pheno_class_${d}.csv" })
}
if (csvFiles.isEmpty()) { println "ERROR: no matching pheno_class_*.csv in ${folder}"; return }

// detections in the CURRENTLY OPEN image, keyed by UUID string
def detById = [:]
getDetectionObjects().each { detById[it.getID().toString()] = it }
println "Open image has ${detById.size()} detections."

int matched = 0, missing = 0
csvFiles.each { file ->
    def reader = new BufferedReader(new FileReader(file))
    reader.readLine()   // header: object_id,broad_lineage,image
    String line
    while ((line = reader.readLine()) != null) {
        def f = line.split(',')
        def obj = detById[f[0]]
        if (obj == null) { missing++; continue }   // belongs to another image — skipped
        obj.setPathClass(getPathClass(f[1]))
        matched++
    }
    reader.close()
}
fireHierarchyUpdate()
println "Set broad lineage on ${matched} detections (${missing} rows for other images / unmatched)."
