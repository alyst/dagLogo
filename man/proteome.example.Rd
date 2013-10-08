\name{proteome.example}
\alias{proteome.example}
\docType{data}
\title{
the subset proteome of fruit fly
}
\description{
the subset proteome of fruit fly
}
\usage{data(proteome.example)}
\format{
  An object of Proteome for fly subset proteome. 
  The format is: A list with one data frame and an character.
  \describe{
    \item{\code{proteome}}{'data.frame':    1406 obs. of  4 variables}
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
    taxId(UniProt.ws) <- 7227
    proteome <- prepareProteome(UniProt.ws)
    proteome@proteome <- proteome@proteome[sample(1:19902, 1406), ]
}
\examples{
data(proteome.example)
head(proteome.example@proteome)
proteome.example@type
}
\keyword{datasets}
