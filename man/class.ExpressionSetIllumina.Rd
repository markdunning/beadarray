\name{ExpressionSetIllumina}
\docType{class}
\alias{class:ExpressionSetIllumina}
\alias{ExpressionSetIllumina-class}


\alias{exprs,ExpressionSetIllumina-method}
\alias{exprs<-,ExpressionSetIllumina,matrix-method}

\alias{se.exprs,ExpressionSetIllumina-method}
\alias{se.exprs<-,ExpressionSetIllumina,matrix-method}

\alias{NoBeads}
\alias{NoBeads,ExpressionSetIllumina-method}
\alias{NoBeads<-}
\alias{NoBeads<-,ExpressionSetIllumina,matrix-method}

\alias{Detection}
\alias{Detection,ExpressionSetIllumina-method}
\alias{Detection<-}
\alias{Detection<-,ExpressionSetIllumina,matrix-method}

\alias{QCInfo}
\alias{QCInfo<-}
\alias{QCInfo<-,ExpressionSetIllumina,list-method}
\alias{QCInfo,ExpressionSetIllumina-method}

\alias{[,ExpressionSetIllumina-method}

\alias{initialize,ExpressionSetIllumina-method}

\alias{combine,ExpressionSetIllumina,ExpressionSetIllumina-method}
\alias{combine,ExpressionSetIllumina,ANY-method}
\alias{show,ExpressionSetIllumina-method}

\alias{getVariance,ExpressionSetIllumina-method}


\title{
  Class to Contain Objects Describing High-Throughput Illumina
  Expression BeadArrays.
}
\description{
  Container for high-throughput assays and experimental
  metadata. \code{ExpressionSetIllumina} class is derived from
  \code{\link[Biobase:class.eSet]{eSet}}, and requires matrices \code{exprs},
  \code{se.exprs}, \code{NoBeads}, \code{Detection} as assay data members.
 }
\section{Extends}{
   Directly extends class \code{\link[Biobase:class.eSet]{eSet}}.
}
\section{Creating Objects}{
  \code{new('ExpressionSetIllumina',
    phenoData = [AnnotatedDataFrame],
    exprs = [matrix],
    se.exprs = [matrix],
    NoBeads = [matrix],
    Detection = [matrix],
    annotation = [character],
    featureData = [AnnotatedDataFrame],
    experimentData = [MIAME],
    ...)
  }
  \code{ExpressionSetIllumina} instances are usually created through
  \code{new("ExpressionSetIllumina", ...)}. Arguments to \code{new}
  include \code{exprs}, \code{se.exprs},
  \code{NoBeads}, \code{Detection}, \code{phenoData},
  \code{experimentData}, and \code{annotation}. 
  \code{phenoData}, \code{experimentData}, and \code{annotation} can be
  missing, in which case they are assigned default values.
}
\section{Slots}{
  Inherited from \code{\link[Biobase:class.eSet]{eSet}}:
  \describe{
    \item{\code{assayData}:}{Contains matrices with equal dimensions, and with 
      column number equal to \code{nrow(phenoData)}. \code{assayData} must 
      contain a matrix \code{exprs} with rows representing features (e.g., genes)
      and columns representing samples, a matrix \code{se.exprs} 
      describing the standard error of each gene, and matrices \code{NoBeads} and \code{Detection} 
      to describe the number of beads used to produce the summary and a
      probability of a gene being expressed above background. The contents of these matrices 
      are	not enforced by the class. Additional matrices of identical size may 
      also be included in \code{assayData}. Class:\code{\link[Biobase:AssayData]{AssayData}}}
    \item{\code{phenoData}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
    \item{\code{experimentData}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
    \item{\code{annotation}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
    \item{\code{protocolData}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
    \item{\code{featureData}:}{annotation for SNPs, usually will contain a 
      \code{CHR} and a \code{MapInfo} column for genomic localization }
  }
 }
\section{Methods}{

  Class-specific methods:
  \describe{
     \item{\code{exprs(ExpressionSetIllumina)}, \code{exprs(ExpressionSetIllumina,matrix)<-}}{Access and
       set elements named \code{exprs} in the \code{AssayData}
       slot.}

     \item{\code{se.exprs(ExpressionSetIllumina)}, \code{se.exprs(ExpressionSetIllumina,matrix)<-}}{Access and
       set elements named \code{se.exprs} in the \code{AssayData}
       slot.}
     \item{\code{NoBeads(ExpressionSetIllumina)}}{Access elements named \code{NoBeads} in the \code{AssayData}
       slot.}
     \item{\code{Detection(ExpressionSetIllumina)}}{Access elements named \code{Detection} in the \code{AssayData}
       slot.}
     \item{\code{getVariance(ExpressionSetIllumina)}}{Calculate
	  bead-type specific variance using \code{se.exprs} and
	  \code{NoBeads} from the \code{AssayData} slot.}
     \item{\code{QCInfo(ExpressionSetIllumina)},\code{QCInfo(ExpressionSetIllumina,list)<-}}{Access elements named \code{QC} in the \code{AssayData}
       slot.}
     \item{\code{object[(index)}:}{Conducts subsetting of matrices and
       phenoData and reporterInfo components }
     \item{\code{combine(ExpressionSetIllumina,ExpressionSetIllumina)}:}{performs \code{union}-like 
		   combination in both dimensions of ExpressionSetIllumina objects}
		 \item{\code{show(ExpressionSetSetIllumina)}}{See
		   \code{\link[Biobase:eSet]{eSet}}}
  }
  
   Derived from \code{\link[Biobase:class.eSet]{eSet}}:
   \describe{
     \item{\code{sampleNames(ExpressionSetSetIllumina)} and \code{sampleNames(ExpressionSetSetIllumina)<-}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
     \item{\code{featureNames(ExpressionSetSetIllumina)}, \code{featureNames(ExpressionSetSetIllumina, value)<-}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
     \item{\code{dims(ExpressionSetSetIllumina)}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
     \item{\code{phenoData(ExpressionSetSetIllumina)}, \code{phenoData(ExpressionSetSetIllumina,value)<-}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
     \item{\code{varLabels(ExpressionSetSetIllumina)}, \code{varLabels(ExpressionSetSetIllumina, value)<-}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
     \item{\code{varMetadata(ExpressionSetSetIllumina)}, \code{varMetadata(ExpressionSetSetIllumina,value)<-}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
     \item{\code{pData(ExpressionSetSetIllumina)}, \code{pData(ExpressionSetSetIllumina,value)<-}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
     \item{\code{varMetadata(ExpressionSetSetIllumina)}, \code{varMetadata(ExpressionSetSetIllumina,value)}}{See \code{\link[Biobase:class.eSet]{eSet}}}
     \item{\code{experimentData(ExpressionSetSetIllumina)},\code{experimentData(ExpressionSetSetIllumina,value)<-}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
%     \item{\code{pubMedIds(ExpressionSetSetIllumina)}, \code{pubMedIds(ExpressionSetSetIllumina,value)}}{See \code{\link[Biobase:class.eSet]{eSet}}}
%     \item{\code{abstract(ExpressionSetSetIllumina)}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
     \item{\code{annotation(ExpressionSetSetIllumina)}, \code{annotation(ExpressionSetSetIllumina,value)<-}}{See \code{\link[Biobase:class.eSet]{eSet}}}
     \item{\code{storageMode(eSet)}, \code{storageMode(eSet,character)<-}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
  }

  Standard generic methods:
  \describe{
    \item{\code{initialize(ExpressionSetSetIllumina)}:}{Object instantiation, used
      by \code{new}; not to be called directly by the user.}
    \item{\code{validObject(ExpressionSetSetIllumina)}:}{Validity-checking method, ensuring
      that \code{call}, \code{callProbability}, \code{G}, and \code{R} are members of
      \code{assayData}. \code{checkValidity(ExpressionSetSetIllumina)} imposes this
      validity check, and the validity checks of \code{Biobase:class.eSet}.}
    \item{\code{show(ExpressionSetSetIllumina)}}{See \code{\link[Biobase:class.eSet]{eSet}}}
    \item{\code{dim(ExpressionSetSetIllumina)}, \code{ncol}}{See \code{\link[Biobase:class.eSet]{eSet}}}
    \item{\code{ExpressionSetSetIllumina[(index)}:}{See \code{\link[Biobase:class.eSet]{eSet}}}
    \item{\code{ExpressionSetSetIllumina$}, \code{ExpressionSetSetIllumina$<-}}{See \code{\link[Biobase:class.eSet]{eSet}}}
  }
}
\author{Mark Dunning, based on Biobase eSet class}
\seealso{\code{\link[Biobase:class.eSet]{eSet}}}
\keyword{classes}
