\name{testDAUresults-class}
\docType{class}
\alias{testDAUresults-class}

\title{Class \code{"testDAUresults"}}
\description{
  An object of class \code{"testDAUresults"} represents background model.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("dagBackground", group="character",
                        difference="matrix",
                        zscore="matrix",
                        pvalue="matrix",
                        background="matrix",
                        motif="matrix",
                        upstream="numeric",
                        downstream="numeric")}.
}
\section{Slots}{
  \describe{
    \item{\code{group}}{Object of class \code{"character"} could be "null", "classic", "charge", "chemistry", "hydrophobicity"}
    \item{\code{difference}}{code{"matrix"} the difference of inputs from background for each amino acid in each position}
    \item{\code{zscore}}{code{"matrix"} z score for each amino acid in each position}
    \item{\code{pvalue}}{code{"matrix"} pvalue for each amino acid in each position}
    \item{\code{background}}{code{"matrix"} background frequencies for each amino acid in each position}
    \item{\code{motif}}{code{"matrix"} inputs frequencies for each amino acid in each position}
    \item{\code{upstream}}{\code{"numeric"} The upstream offset from anchor}
    \item{\code{downstream}}{\code{"numeric"} The downstream offset from anchor}
  }
}
\keyword{classes}
