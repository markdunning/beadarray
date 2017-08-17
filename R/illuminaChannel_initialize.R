

setMethod("initialize", "illuminaChannel",
          function(.Object,  transFun, outlierFun, exprFun, varFun,channelName) {
               .Object@name = channelName
                .Object@transFun = list(transFun)
               .Object@outlierFun = list(outlierFun)
               .Object@exprFun = list(exprFun)
		.Object@varFun = list(varFun)
.Object})

