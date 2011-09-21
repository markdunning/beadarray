\name{plotByProbeFactor}
\alias{plotByProbeFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a particular probe factor
}
\description{
Produces a boxplot of expression level against a probe-wise factor of interest.
}
\usage{
plotByProbeFactor(data, probeFactor = "PROBEQUALITY")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An ExpressionSetIllumina object
}
  \item{probeFactor}{
A column in the featureData slot for the ExpressionSetIllumina object, or vector of equivalent length
}
}
\details{

The column in the featureData is identified and then extracted as vector whose length is the same as the number of probes. Alternatively, this vector could be passed directly to the function.

The vector is then merged with the expression matrix of the ExpressionSetIllumina and reshaped into a form that is convenient for ggplot graphics system. A boxplot is then produced with expression level of a probe as a response and the specified vector as a factor. 


}


\value{
A ggplot object

}
\author{
Mark Dunning
}
