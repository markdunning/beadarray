## Class for storing bead-level data
setClass("beadLevelData",
          representation(beadData = "list",
            sectionData = "list",
            experimentData = "list",
            history = "character"))
            
## The previous class for storing bead-level data, now deprecated.           
## Maintained so previously saved objects can be converted to the new structure
setClass("BeadLevelList",
          representation(beadData = "list",
            arrayInfo = "list",
            phenoData = "AnnotatedDataFrame",
            annotation="character"))

setClass("illuminaChannel",
         representation(transFun ="list",outlierFun="list", exprFun="list", varFun = "list",name="character"),
        
)

setClass("ExpressionSetIllumina",
         representation(QC = "AnnotatedDataFrame", channelData="list"),
         contains="eSet"
)

setClass("limmaResults",
         representation(DesignMatrix = "matrix", ArrayWeights = "vector", ContrastMatrix = "matrix", annotation="character"),
         contains="eSet"
)

setClass("beadRegistrationData",
          representation(
                        layout = "list",                                                
                        registrationData = "list",
                        coordinateData = "list",
                        cornerData = "list",
                        p95 = "numeric",
                        imageLocations = "character", 
                        metrics = "data.frame")
)  
