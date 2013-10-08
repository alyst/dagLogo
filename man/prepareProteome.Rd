\name{prepareProteome}
\alias{prepareProteome}
\title{
prepare proteome for background building
}
\description{
prepare proteome from UniProt webserver or a fasta file
}
\usage{
prepareProteome(UniProt.ws, fasta, species="unknown")
}
\arguments{
  \item{UniProt.ws}{an object of UniProt.ws}
  \item{fasta}{fasta file name or an object of AAStringSet}
  \item{species}{an character to assign the species of the proteome}
}
\value{
an object of Proteome which contain protein sequence information
}
\author{
Jianhong Ou
}
\seealso{
\code{\link{formatSequence}}, \code{\link{buildBackgroundModel}}
}
\examples{
    if(interactive()){
        taxId(UniProt.ws) <- 7227
        proteome <- prepareProteome(UniProt.ws, species="Drosophila melanogaster")
    }
}

\keyword{ misc }