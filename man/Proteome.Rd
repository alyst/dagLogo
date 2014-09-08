\name{Proteome-class}
\docType{class}
\alias{Proteome-class}

\title{Class \code{"Proteome"}}
\description{
  An object of class \code{"Proteome"} represents proteome of a given species.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form 
  \code{new("Proteome", proteome, type, species)}.
}
\section{Slots}{
  \describe{
    \item{\code{proteome}}{Object of class \code{"data.frame"} the proteome of a given species, 
    should include ids and peptide sequences.}
    \item{\code{type}}{code{"character"} indicates how the object is prepared, 
    could be "fasta" or "UniProt"}
    \item{\code{species}}{\code{"character"} the species}
  }
}
\keyword{classes}
