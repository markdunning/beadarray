\name{readSampleSheet}
\alias{readSampleSheet}
\alias{sampleSheet}
\alias{sampleSheet,ExpressionSetIllumina-method}
\alias{sampleSheet,beadLevelData-method}
\alias{sampleSheet<-}
\alias{sampleSheet<-,ExpressionSetIllumina,data.frame-method}
\alias{sampleSheet<-,beadLevelData,data.frame-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read a Sample sheet for a BeadArray experiment
}
\description{
Read a comma-separated sample sheet into beadarray, which can then be used to annotate the data that beadarray reads in and add relevant metadata to the analysis.
}
\usage{
readSampleSheet(sheet = "sampleSheet.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sheet}{
A file name for the sample sheet
}
}
\details{
The sample sheet should describe all sections are expected to be part of the experiment, along with any relevant metadata. We assume the GenomeStudio sample sheet specification, with [Data] and [Header] sections. The Data section is where the sample metadata are defined. One row is required for each section, and the columns Sentrix\_ID, Sentrix\_Position, Sample\_Name, Sample\_Group must be present. The Sentrix\_ID and Sentrix\_Position columns give the name of the chip (usually a 10 digit number) and slot (a letter A to L) where a given sample was hybridised to. The Sample\_Group and Sample\_Name columns can take any valid R name. If there are muliple images per section, for example for Mouse WG-6 data, which have bead-level data suffixed by A\_1, A\_2, B\_1, B\_2 etc, only one entry is required in the sample sheet for sections A and B. beadarray will duplicate the metadata as neccessary.

}
\value{
a list comprised of a tabular sample table and any header information that is present.

}
\author{
Mark Dunning
}
