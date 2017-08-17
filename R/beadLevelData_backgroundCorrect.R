
backgroundCorrectSingleSection = function(BLData, array=1 ,fg="Grn" , bg="GrnB" ,newName = "Grn.bc"){

##Check that the value of array is valid

	fgI = getBeadData(BLData, array=array, what=fg)
	bgI = getBeadData(BLData, array=array, what=bg)	

	data = fgI - bgI

	BLData = insertBeadData(BLData, array=array, data = data, what = newName)

	

	BLData
}


