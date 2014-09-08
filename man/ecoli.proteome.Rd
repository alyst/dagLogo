\name{ecoli.proteome}
\alias{ecoli.proteome}
\docType{data}
\title{
the subset proteome of Escherichia coli
}
\description{
the subset proteome of Escherichia coli
}
\usage{data(ecoli.proteome)}
\format{
  An object of Proteome for Escherichia coli proteome. 
  The format is: A list with one data frame and an character.
  \describe{
    \item{\code{proteome}}{'data.frame':     obs. of  4 variables}
    \item{\code{type}}{'character':   "UniProt"}
  }
  The format of proteome is
  \describe{
        \item{\code{ENTREZ_GENE}}{a character vector, records entrez gene id}
        \item{\code{SEQUENCE}}{a character vector, peptide sequences}
        \item{\code{ID}}{a character vector, Uniprot ID}
        \item{\code{LEN}}{a character vector, length of peptides}
    }
}
\details{
used in the examples
  Annotation data obtained by:
    library(UniProt.ws)
    taxId(UniProt.ws) <- 562
    proteome <- prepareProteome(UniProt.ws, species="Escherichia coli")
}
\examples{
data(ecoli.proteome)
head(ecoli.proteome@proteome)
ecoli.proteome@type
}
\keyword{datasets}
