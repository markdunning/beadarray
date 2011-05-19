\name{show-method}
\docType{methods}
\alias{show,beadLevelData-method}
\alias{show,ExpressionSetIllumina-method}
\alias{show,beadRegistrationData-method}
\title{Display object summary}
\description{
	Prints a summary of an objects contents.
}
\usage{

\S4method{show}{beadLevelData}(object)

\S4method{show}{ExpressionSetIllumina}(object)

}
\arguments{
\item{object}{An object of class \code{\link[=beadLevelData-class]{beadLevelData}} or \code{\link[=ExpressionSetIllumina-class]{ExpressionSetIllumina}} }
}

\details{
    \code{show} is commonly invoked by simply entering the name of an object.  Calling it on the classes defined in beadarray will print a summary of the object contents, with the actual output dependent on the class of the object.
    
    Output for the \code{\link[=beadLevelData-class]{beadLevelData}} class is broken down into three sections: experiment information, data that relate to each array section and data for individual beads.  The full information relating to the first two groups will be display, with only a short summary of the per-bead information shown (currently 5 beads from the first section).
    
    The \code{\link[=ExpressionSetIllumina-class]{ExpressionSetIllumina}} class is based up the eSet class and the output from \code{show} is closely related, with a short summary of the contents of available slots.
}
    

\author{Mark Dunning}


\keyword{methods}
