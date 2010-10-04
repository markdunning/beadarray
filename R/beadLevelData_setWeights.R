setWeights = function(BLData, wts, array, combine = FALSE, wtName = "wts"){

	if (length(array) == 1)
	{
		##one array
		if(combine)
		{
			#check for existing weights
			if(!wtName %in% colnames(BLData[[array]]))
			{
				BLData = insertBeadData(BLData, array=array, what=wtName, data = wts[[array]])
			}
			else
			{
				wtCol = grep(wtName, colnames(BLData[[array]]))

				BLData = insertBeadData(BLData, array=array, what=wtName,data = pmin(wts[[array]],BLData[[array]][,wtCol]))
			}
		}
		else
		{
			BLData = insertBeadData(BLData, array=array, what=wtName, data = wts[[array]])
		}
	}
	else
	{
		##multiple arrays
		if(combine)
		{
			for(i in array)
			{
				#check for existing weights
				if(!wtName %in% colnames(BLData[[array]]))
				{BLData = insertBeadData(BLData, array=i, what=wtName, data = wts[[i]])}
				else{
					wtCol = grep(wtName, colnames(BLData[[array]]))

					BLData = insertBeadData(BLData, array=i, what=wtName, data = pmin(wts[[i]],BLData[[i]][,wtCol]))
				}
			}		
		}
		else
		{
			for(i in array)
			{
				BLData = insertBeadData(BLData, array=i, what=wtName, data = wts[[i]])
			}
		}
	}

BLData

}
