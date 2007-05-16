
setClass("ExpressionSetIllumina",
         representation(QC ="list"),
         contains="eSet"
)

# Define new class, modified from the 'cytoSet' class in prada
setClass("BeadLevelList",
          representation(beadData="environment",
             phenoData="AnnotatedDataFrame",
             arrayInfo ="list",
#             probeindex="matrix",
             annotation="character",
             beadAnno="data.frame",
             scanMetrics="data.frame"))
