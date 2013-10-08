\name{testDAU}
\alias{testDAU}
\title{
DAU test
}
\description{
Performs DAU test
}
\usage{
testDAU(dagPeptides, dagBackground, 
                    group=c("null", "classic", "charge", "chemistry", "hydrophobicity"))
}
\arguments{
  \item{dagPeptides}{an object of dagPeptides, output of \code{\link{fetchSequence}} or f\code{\link{formatSequence}}}
  \item{dagBackground}{an object of dagBackground, output of \code{\link{buildBackgroundModel}}}
  \item{group}{could be "null", "classic", "charge", "chemistry", "hydrophobicity"}
}
\value{
an object of testDAUresults ready for plotting
}
\author{
Jianhong Ou, Julie Zhu
}
\examples{
    data("seq.example")
    data("proteome.example")
    bg <- buildBackgroundModel(seq.example, proteome=proteome.example)
    t <- testDAU(seq.example, bg)
}

\keyword{ misc }