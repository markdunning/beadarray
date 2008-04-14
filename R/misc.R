beadarrayUsersGuide <- function(view=TRUE, topic="beadlevel")
#       function modified from limmaUsersGuide() in limma package
{
        if(topic!="beadlevel" || topic!="beadsummary") {
                cat("\'topic\' must be one of \"beadlevel\" or \"beadsummary\".  Setting topic=\"beadlevel\"\n") 
                topic="beadlevel"
        }
	f = system.file("doc", paste(topic,".pdf",sep=""), package="beadarray")
        if(view) {
                if(.Platform$OS.type == "windows")
                     shell.exec(f)
                else
                     system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
                }
        return(f)
}