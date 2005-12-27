"readBeadSummaryTargets" <- function (file = "beadSummaryTargets.txt", path = NULL,  header=T,sep="" )
{
if (!is.null(path)) file = file.path(path, file)

read.table(file, header=header, sep=sep)

}