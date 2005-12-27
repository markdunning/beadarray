\name{plotCluster}
\alias{plotCluster}
\alias{probeCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Clustering of bead position }
\description{
Performs a clustering on the x and y positions of beads of the same type and displays result
}
\usage{
plotCluster(coords, main = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coords}{vector of numeric values specifying x and y positions of the beads that we're interested in (this can be found using the 
	'getProbeCoords' function}
\item{main}{specify a title for the plot}

}
\details{
Co-ordinates are used to create a distance matrix via the \code{dist} command and then clustered using the ward method. After the clustering we 
choose to divide the co-ordinates in 8 clusters. 

}
\value{
A dendrogram is plotted
}

\author{Mark Dunning}
\keyword{hplot}