setClassUnion("AssayDataOrNull", c("AssayData", "NULL"))
setClassUnion("matrixOrNull", c("matrix", "NULL"))

setClass("ExpressionSetIllumina",
         representation(QC = "AssayDataOrNull"),
         contains="eSet"
)

# Define new class, modified from the 'cytoSet' class in prada
setClass("BeadLevelList",
          representation(beadData="environment",
             phenoData="AnnotatedDataFrame",
             arrayInfo ="list",
#             probeindex="matrix",
             beadAnno="data.frame",
             scanMetrics="data.frame"))
