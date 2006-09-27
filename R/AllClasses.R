setClassUnion("AssayDataOrNull", c("AssayData", "NULL"))


setClass("ExpressionSetIllumina",
         representation(QC = "AssayDataOrNull"),
         contains="eSet"

)

setClass("BeadLevelList",

	representation(
	G="matrix",
	Gb="matrix",
	GrnY ="matrix",
	GrnX="matrix",
	ProbeID="matrix",
	targets ="data.frame"
)

)

