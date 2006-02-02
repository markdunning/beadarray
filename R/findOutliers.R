"findOutliers" <-
function(BLData, probe, array, log=FALSE, n=3){


inten = getProbeIntensities(BLData, probe, array, log)
probe_ids = getProbeIndices(BLData, probe, array)

#nas will be a list of beads which have NA intensity
nas=NULL

if(length(is.na(inten))>0){

nas = probe_ids[is.na(inten) | inten < 0]
probe_ids = probe_ids[!is.na(inten)]
inten = inten[!is.na(inten)]


}




if(log){

raw_inten = getProbeIntensities(BLData, probe,array, log=FALSE)
raw_inten = raw_inten[!is.na(raw_inten)]

}

else{
raw_inten = inten

}



outliers  = NULL

m = mean(inten, na.rm=TRUE, trim=0.5)
ma = mad(inten, na.rm=TRUE)

outliers=probe_ids[inten > m + n *ma | inten < m - n*ma | raw_inten < 0]

return = list()

return$outliers = append(outliers, nas)

return$valid = probe_ids[inten < m + n*ma & inten > m - n*ma & raw_inten > 0]

return

}
