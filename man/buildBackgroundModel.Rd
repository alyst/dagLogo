\name{buildBackgroundModel}
\alias{buildBackgroundModel}
\title{
build background model
}
\description{
build background model for dag test
}
\usage{
buildBackgroundModel(dagPeptides, 
                                 bg=c("wholeGenome", "inputSet", "nonInputSet"),
                                 model=c("any", "anchored"),
                                 targetPosition=c("any", "Nterminus", "Cterminus"),
                                 uniqueSeq=TRUE,
                                 permutationSize=30L,
                                 rand.seed=1,
                                 replacement=FALSE,
                                 proteome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dagPeptides}{an object of dagPeptides, output of \code{\link{fetchSequence}} or \code{\link{formatSequence}}}
  \item{bg}{could be "wholeGenome", "inputSet" or "nonInputSet"}
  \item{model}{could be "any" or "anchored"}
  \item{targetPosition}{could be "any", "Nterminus" or "Cterminus"}
  \item{uniqueSeq}{should the background sequence be unique?}
  \item{permutationSize}{how many times should it samples}
  \item{rand.seed}{random seed}
  \item{replacement}{Should sampling be with replacement?}
  \item{proteome}{an object of Proteome, output of \code{\link{prepareProteome}}}
}
\details{
    The background could be generated from wholeGenome, inputSet or nonInputSet.
    whole genome: randomly select subsequences from the whole genome with each subsequence 
    containing amino acids with same width of input sequences.
    anchored whole genome: randomly select subsequences from the whole genome with 
    each subsequence containing amino acids with same width of input sequences
    where the middle amino acids must contain anchor amino acid, e.g., K, which is
    specified by user.
    input set: same to whole genome, but only use protein sequence from input id and not 
    including the site specified in input sequences
    anchored input set: same to anchored whole genome, but only use protein sequences from input id, 
    and not including the site specified in input sequences.
    non-input set: whole genome - input set.
    anchored non-input set: whole genome - input set and the middle amino acids must contain anchor amino acid.
}
\value{
    an object of dagBackground which contains background and permutationSize.
}
\author{
    Jianhong Ou, Alexey Stukalov, Julie Zhu
}
\seealso{
    \code{\link{prepareProteome}}
}
\examples{
    data("seq.example")
    data("proteome.example")
    bg <- buildBackgroundModel(seq.example, proteome=proteome.example)
}
\keyword{ misc }