# 
# 
# 
# 
# applySummary = function(x, outlierFunction, summaryFunctions){
# 
# res = NULL
# 
# i=1
# 
# x = outlierFunction(x)
# 
# 	for(fn in summaryFunctions){
# 
# 		res[i] = fn(x)
# 		i = i+1
# 
# 	}
# 
# res
# 
# }
# 
# setClass("illuminaChannel",
#          representation(transFun ="list",outlierFun="list", summaryFun="list", name="character"),
#         
# )
# 
# 
# setMethod("initialize", "illuminaChannel",
#           function(.Object,  transFun, outlierFun, summaryFun, channelName) {
#                .Object@name = channelName
#                 .Object@transFun = list(transFun)
#                .Object@outlierFun = list(outlierFun)
#                .Object@summaryFun = list(summaryFun)
# 			
#                .Object})
# 
# createBeadSummaryData = function(BLData, channelList){
# 
# nArrays = length(arrayNames(BLData))
# 
# ####Asign an NChannelSet with the correct number of channels, rows, columns etc....
# 
# for(i in 1:nArrays){
# 
# cat("Processing Array", i, "\n")
# 
# for(j in 1:length(channelList)){
# 
# 
# transformFun = channelList[[j]]@transFun[[1]]
# outlierMethod = channelList[[j]]@outlierFun[[1]]
# summaryFunctions = channelList[[j]]@summaryFun[[1]]
# 
# ####Eventually this should add data to the NChannelSet....
# bsd <- t(sapply(transformFun(BLData,i), function(x) applySummary(x, outlierMethod, summaryFunctions)))
# 
# }
# 
# }
# 
# bsd
# 
# }
# 
# #####################An example script###########################################
# 
# 
# ###Read your favourite bead-level data, or use data(BLData)
# 
# ###BLData = readIllumina(useImages=FALSE)
# 
# 
# ####Define the functions
# 
# summaryFunctions=c(mean, sd, length)
# 
# ####Define an outlier method
# 
# illuminaOutlierMethod = function(x){
# 
# M = median(x)
# 
# md = mad(x)
# 
# x[x < M+3*md & x > M-3*md]
# 
# }
# 
# 
# ###A transformation funtion that gives a list of quantities to be summarized indexed by ProbeID
# 
# singleChannelTransform = function(BLData,array){
# 
# ###You also need code to remove masked beads etc
# 
# pids = getArrayData(BLData, what="ProbeID", array=array)
# inten = getArrayData(BLData, what="G", array=array)
# split(inten, pids)
# 
# }
# 
# 
# ###A transformation to get log-ratios might look like this
# 
# logRatioTransform = function(BLData){
# pids = getArrayData(BLData, what="ProbeID", array=1)
# Ginten = getArrayData(BLData, what="G", array=1)
# Rinten = getArrayData(BLData, what="G", array=1)
# 
# lr = R - G
# 
# split(lr, pids)
# 
# }
# 
# 
# 
# 
# singleChannel  = new("illuminaChannel", singleChannelTransform, illuminaOutlierMethod, summaryFunctions, "G")
# 
# system.time(bsd <- createBeadSummaryData(BLData, list(singleChannel)))
# 
# 
# ####For two-colour data, we may want to do the following
# 
# logRatio = new("illuminaChannel", logRatioTransform, illuminaOutlierMethod, summaryFunctions, "M")
# 
# redChannel = new("illuminaChannel", singleChannelTransform, illuminaOutlierMethod, summaryFunctions, "R")
# greenChannel = new("illuminaChannel", singleChannelTransform, illuminaOutlierMethod, summaryFunctions, "G")
# 
# ##bsd2 = createBeadSummaryData(BLData, list(singleChannel, redChannel, logRatio))
# 
# 
# 
# 
