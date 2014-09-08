\name{dagPeptides-class}
\docType{class}
\alias{dagPeptides-class}

\title{Class \code{"dagPeptides"}}
\description{
  An object of class \code{"dagPeptides"} represents the information of peptides.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("dagPeptides", data, peptides, upstreamOffset, downstreamOffset, type)}.
}
\section{Slots}{
  \describe{
    \item{\code{data}}{Object of class \code{"data.frame"} The details of the input sequences. 
    It includes the columns: IDs, anchorAA (anchor Amino Acid), anchorPos (anchor Position),
    peptide (protein peptide), anchor, upstream, downstream (peptides in given upstream and downstream offset from anchor)}
    \item{\code{peptides}}{code{"matrix"} The input peptides. Each column contains one peptide in that position}
    \item{\code{upstreamOffset}}{\code{"numeric"} The upstream offset from anchor}
    \item{\code{downstreamOffset}}{\code{"numeric"} The downstream offset from anchor}
    \item{\code{type}}{\code{"charactger"} ID type of inputs}
  }
}
\keyword{classes}
