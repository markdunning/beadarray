\name{makeQCProfile}
\alias{makeControlProfile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve control beads
}

\description{
Function to retrieve the bead IDs for a particular annotation platform
}

\usage{

makeControlProfile(annoName, excludeERCC = TRUE)

}



%- maybe also 'usage' for other objects documented here.
\arguments{

\item{annoName}{An abberviated name for the annotation packaget. e.g. Humanv3, MouseV2}
\item{excludeERCC}{Ignore any controls with REPORTERGROUPNAME beginning with ERCC. These are a large collection of spike controls on the Humanv4 chip and can make some control tables too big if they are included.}
}

\details{

The function will use the relevant annotation package (if it is installed) and return all IDs that are present in the REPORTERGROUPNAME mapping. Note that this mapping is only available in versions  1.12.0 and greater.


}

\value{
A matrix with one row per control bead-type giving the bead-level ID (ArrayAddressID) and control type.
}

\examples{

makeControlProfile("Humanv3")

}

\author{
Mark Dunning
}


