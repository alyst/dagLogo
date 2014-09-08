\name{formatSequence}
\alias{formatSequence}
\title{
    prepare an object of dagPeptides from sequences
}
\description{
    prepare an object of dagPeptides from sequences
}
\usage{
formatSequence(seq, proteome, upstreamOffset, downstreamOffset)
}
\arguments{
  \item{seq}{a vector of character, amino acid sequences}
  \item{proteome}{an object of Proteome, output of \code{\link{prepareProteome}}}
  \item{upstreamOffset}{an integer, upstream offset position}
  \item{downstreamOffset}{an integer, downstream offset position}
}
\value{
return an object of dagPeptides, which is a list contains:
    data, peptides, upstreamOffset, downstreamOffset and type information
}\author{
    Jianhong Ou, Julie Zhu
}

\seealso{
\code{\link{fetchSequence}}
}
\examples{
if(interactive()){
    dat <- unlist(read.delim(system.file("extdata", 
                            "grB.txt", package="dagLogo"), 
                            header=F, as.is=TRUE))
    proteome <- prepareProteome(fasta=system.file("extdata", 
                                              "HUMAN.fasta",
                                              package="dagLogo"))
    seq <- formatSequence(dat, proteome)
}
}
\keyword{ misc }