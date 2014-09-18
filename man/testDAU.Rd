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
                    group=c("null", "classic", "charge", "chemistry", "hydrophobicity"),
                    bgNoise=NA)
}
\arguments{
  \item{dagPeptides}{an object of dagPeptides, output of \code{\link{fetchSequence}} or f\code{\link{formatSequence}}}
  \item{dagBackground}{an object of dagBackground, output of \code{\link{buildBackgroundModel}}}
  \item{group}{could be "null", "classic", "charge", "chemistry", "hydrophobicity"}
  \item{bgNoise}{if it is not NA, test will using a background by Dirichlet(1)-distributed random frequencies with weight bg.noise. 
  The value of bgNoise should be a number in the range of 0 to 1, eg. 0.05}
}
\value{
an object of testDAUresults ready for plotting
}
\author{
Jianhong Ou, Alexey Stukalov, Julie Zhu
}
\examples{
    data("seq.example")
    data("proteome.example")
    bg <- buildBackgroundModel(seq.example, proteome=proteome.example)
    t <- testDAU(seq.example, bg, bgNoise=0.05)
}

\keyword{ misc }