 
setClass("ExpressionSetIllumina",     
         
   contains = "eSet"
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

