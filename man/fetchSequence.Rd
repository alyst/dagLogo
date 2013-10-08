\name{fetchSequence}
\alias{fetchSequence}
\title{
fetch sequence by id
}
\description{
fetch amino acid sequence by given identifiers via biomaRt or proteome prepared by \code{\link{prepareProteome}}
}
\usage{
fetchSequence(IDs, type="entrezgene", anchorAA=NULL, anchorPos,
                            mart, proteome, upstreamOffset, downstreamOffset)
}
\arguments{
  \item{IDs}{A vector of Identifiers to retrieve peptides}
  \item{type}{type of identifiers}
  \item{anchorAA}{a vector of character, anchor Amino Acid}
  \item{anchorPos}{a vector of character or numeric, anchor position, for example, K121}
  \item{mart}{an object of Mart}
  \item{proteome}{an object of Proteome, output of \code{\link{prepareProteome}}}
  \item{upstreamOffset}{an integer, upstream offset position}
  \item{downstreamOffset}{an integer, downstream offset position}
}
\value{
    return an object of \code{\linkS4class{dagPeptides}}
}
\author{
    Jianhong Ou, Julie Zhu
}

\seealso{
\code{\link{formatSequence}}
}
\examples{
    mart <- useMart("ensembl", "dmelanogaster_gene_ensembl")
    dat <- read.csv(system.file("extdata", "dagLogoTestData.csv", package="dagLogo"))
    seq <- fetchSequence(as.character(dat$entrez_geneid[1:5]), 
                anchorPos=as.character(dat$NCBI_site[1:5]), 
                mart=mart, 
                upstreamOffset=7, 
                downstreamOffset=7)
}
\keyword{ misc }