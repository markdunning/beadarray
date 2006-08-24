 
setClass("ExpressionSetIllumina",     
         
   contains = "eSet"
)

setClass("ExpressionQCSet", representation(
         Signal = "matrix",
         Var = "matrix",
         Detection = "matrix"

                )
         )
