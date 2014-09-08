\name{seq.example}
\alias{seq.example}
\docType{data}
\title{
example object of dagPeptides
}
\description{
example object of dagPeptides
}
\usage{data(seq.example)}
\format{
   An object of dagPeptides. 
  The format is: A list.
  \describe{
    \item{\code{data}}{'data.frame':    732 obs. of  7 variables}
    \item{\code{peptides}}{'matrix':   amnio acid in each position}
    \item{\code{upstreamOffset}}{an integer, upstream offset position}
    \item{\code{downstreamOffset}}{an integer, downstream offset position}
    \item{\code{type}}{"character", type of identifiers}
  }
  The format of data is
  \describe{
        \item{\code{IDs}}{a character vector, input identifiers}
        \item{\code{anchorAA}}{a character vector, anchor amino acid provided in inputs}
        \item{\code{anchorPos}}{a numeric vector, anchor position in the protein}
        \item{\code{peptide}}{a character vector, peptide sequences}
        \item{\code{anchor}}{a character vector, anchor amino acid in the protein}
        \item{\code{upstream}}{a character vector, upstream peptides}
        \item{\code{downstream}}{a character vector, downstream peptides}
    }
}
\details{
used in the examples
  seq obtained by:
    mart <- useMart("ensembl", "dmelanogaster_gene_ensembl")
    dat <- read.csv(system.file("extdata", "dagLogoTestData.csv", package="dagLogo"))
    seq <- fetchSequence(as.character(dat$entrez_geneid), 
                anchorPos=as.character(dat$NCBI_site), 
                mart=mart, 
                upstreamOffset=7, 
                downstreamOffset=7)
}
\examples{
data(seq.example)
head(seq.example@peptides)
seq.example@upstreamOffset
seq.example@downstreamOffset
}
\keyword{datasets}
