\name{ExpressionSetIllumina}
\docType{class}
\alias{class:ExpressionSetIllumina}
\alias{ExpressionSetIllumina-class}
\alias{exprs<-,ExpressionSetIllumina,matrix-method}
\alias{exprs,ExpressionSetIllumina-method}
\alias{[,ExpressionSetIllumina-method}
\alias{initialize,ExpressionSetIllumina-method}
\alias{combine,ExpressionSetIllumina,ExpressionSetIllumina-method}
\alias{combine,ExpressionSetIllumina,ANY-method}

\title{
  Class to Contain Objects Describing High-Throughput Illumina
  Expression BeadArrays.
}


\description{
  Container for high-throughput assays and experimental
  metadata. \code{ExpressionSetIllumina} class is derived from
  \code{\link[Biobase:eSet]{eSet}}, and requires matrices \code{exprs}, \code{BeadStDev}, \code{NoBeads},
  \code{Detection} as assay data members.  
  
 } 


\section{Extends}{
   Directly extends class \code{\link[Biobase:eSet]{eSet}}.
}
\section{Creating Objects}{

  \code{new('ExpressionSetIllumina',
    phenoData = [AnnotatedDataFrame],
    experimentData = [MIAME],
    annotation = [character],
    exprs = [matrix],
    BeadStDev = [matrix],
    NoBeads = [matrix],
    Detection = [matrix],
    reporterInfo = [data.frameOrNULL],
    ...)
  }

  \code{ExpressionSetIllumina} instances are usually created through
  \code{new("ExpressionSetIllumina", ...)}. Arguments to \code{new}
  include \code{exprs}, \code{BeadStDev},
  \code{NoBeads}, \code{Detection}, \code{phenoData},
  \code{experimentData}, and \code{annotation}. 
  \code{phenoData}, \code{experimentData}, and \code{annotation} can be
  missing, in which case they are assigned default values.
}
\section{Slots}{
  Inherited from \code{Biobase:eSet}:
  \describe{
    \item{\code{assayData}:}{Contains matrices with equal dimensions, and with 
      column number equal to \code{nrow(phenoData)}. \code{assayData} must 
      contain a matrix \code{exprs} with rows representing features (e.g., genes) 
      and columns representing samples, a matrix \code{BeadStDev} 
      describing the standard error of each gene, and matrices \code{NoBeads} and \code{Detection} 
      to describe the number of beads used to produce the summary and a
      probability of a gene being expressed above background. The contents of these matrices 
      are	not enforced by the class. Additional matrices of identical size may 
      also be included in \code{assayData}. Class:\code{\link[Biobase:AssayData]{AssayData}}}
    \item{\code{phenoData}:}{See \code{\link[Biobase:eSet]{eSet}}}
    \item{\code{experimentData}:}{See \code{\link[Biobase:eSet]{eSet}}}
    \item{\code{annotation}:}{See \code{\link[Biobase:eSet]{eSet}}}
    \item{\code{reporterInfo}:}{annotation for SNPs, usually will contain a 
      \code{CHR} and a \code{MapInfo} column for genomic localization }
  }
}
\section{Methods}{

  Class-specific methods:
  \describe{
     \item{\code{se.exprs(ExpressionSetIllumina)}, \code{se.exprs(ExpressionSetIllumina,matrix)<-}}{Access and
       set elements named \code{BeadStDev} in the \code{AssayData}
       slot.}



     
     \item{\code{reporterInfo(ExpressionSetIllumina)}, \code{reporterInfo(ExpressionSetIllumina,data.frame)<-}}{
		   access to the \code{reporterInfo} slot}  
     \item{\code{object[(index)}:}{Conducts subsetting of matrices and
       phenoData and reporterInfo components }
     \item{\code{combine(ExpressionSetIllumina,ExpressionSetIllumina)}:}{performs \code{union}-like 
		   combination in both dimensions of ExpressionSetIllumina objects}
	       }

	       
Derived from \code{\link[Biobase:eSet]{eSet}}:
   \describe{
     \item{\code{sampleNames(ExpressionSetIllumina)} and \code{sampleNames(ExpressionSetIllumina)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{featureNames(ExpressionSetIllumina)}, \code{featureNames(ExpressionSetIllumina, value)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{dims(ExpressionSetIllumina)}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{phenoData(ExpressionSetIllumina)}, \code{phenoData(ExpressionSetIllumina,value)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{varLabels(ExpressionSetIllumina)}, \code{varLabels(ExpressionSetIllumina, value)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{varMetadata(ExpressionSetIllumina)}, \code{varMetadata(ExpressionSetIllumina,value)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{pData(ExpressionSetIllumina)}, \code{pData(ExpressionSetIllumina,value)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{varMetadata(ExpressionSetIllumina)}, \code{varMetadata(ExpressionSetIllumina,value)}}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{experimentData(ExpressionSetIllumina)},\code{experimentData(ExpressionSetIllumina,value)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{pubMedIds(ExpressionSetIllumina)}, \code{pubMedIds(ExpressionSetIllumina,value)}}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{abstract(ExpressionSetIllumina)}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{annotation(ExpressionSetIllumina)}, \code{annotation(ExpressionSetIllumina,value)<-}}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{storageMode(eSet)}, \code{storageMode(eSet,character)<-}:}{See \code{\link[Biobase:eSet]{eSet}}}
     \item{\code{reporterNames(ExpressionSetIllumina)}, \code{reporterNames(ExpressionSetIllumina,value)<-}:}{DEPRECATED}
  }

  Standard generic methods:
  \describe{
    \item{\code{initialize(ExpressionSetIllumina)}:}{Object instantiation, used
      by \code{new}; not to be called directly by the user.}
    \item{\code{validObject(ExpressionSetIllumina)}:}{Validity-checking method, ensuring
      that \code{exprs}, \code{BeadStDev}, \code{NoBeads}, and \code{Dectection} are members of
      \code{assayData}. \code{checkValidity(ExpressionSetIllumina)} imposes this
      validity check, and the validity checks of eSet}
    \item{\code{show(ExpressionSetIllumina)}}{See \code{\link[Biobase:eSet]{eSet}}}
    \item{\code{dim(ExpressionSetIllumina)}, \code{ncol}}{See \code{\link[Biobase:eSet]{eSet}}}
    \item{\code{ExpressionSetIllumina[(index)}:}{See \code{\link[Biobase:eSet]{eSet}}}
    \item{\code{ExpressionSetIllumina$}, \code{ExpressionSetIllumina$<-}}{See \code{\link[Biobase:eSet]{eSet}}}
  }
}

}
}
\author{Mark Dunning, based on Biobase eSet class }

\seealso{
  \code{\link{Biobase:eSet}}}

\keyword{classes}
