readBGX = function(filename, path=".", sep="\t",
       quote="", header=TRUE, probeStart="\\[Probes\\]",
       controlStart="\\[Controls\\]", ...) {
  tmp = readLines(file.path(path, filename))
  skip = grep(probeStart, tmp)
  end = grep(controlStart, tmp)
  if(header)
    nrows = end-skip-2
  else
    nrows = end-skip-1
  probeAnno = read.table(file.path(path, filename), sep=sep, header=header, skip=skip, nrows=nrows, quote=quote, ...)
  controlAnno = read.table(file.path(path, filename), sep=sep, header=header, skip=end, ...)
  tmp = probeAnno[1:nrow(controlAnno),]
  ind = match(colnames(controlAnno), colnames(probeAnno))
  tmp[,-ind[!is.na(ind)]] = ""
  tmp[,ind[!is.na(ind)]] = controlAnno[,!is.na(ind)]
  Status = rep("gene", nrow(probeAnno)+nrow(controlAnno))
  Status[(nrow(probeAnno)+1):(nrow(probeAnno)+nrow(controlAnno))] = as.character(controlAnno$Reporter_Group_Name)
  anno = rbind(probeAnno, tmp)
  anno = cbind(anno,Status)
  anno
}
