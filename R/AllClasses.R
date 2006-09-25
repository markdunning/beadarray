setClassUnion("AssayDataOrNull", c("AssayData", "NULL"))
setClassUnion("matrixOrNull", c("matrix", "NULL"))

setClass("ExpressionSetIllumina",
         representation(QC = "AssayDataOrNull"),
         contains="eSet"

)

setClass("BeadLevelList",

	representation(
	G="matrix",
	Gb="matrix",
	R="matrixOrNull",
	Rb="matrixOrNull",
	GrnY ="matrix",
	GrnX="matrix",
	ProbeID="matrix",
	targets ="data.frame"
)

)

