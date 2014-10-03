\name{dagLogo}
\alias{dagLogo}
\title{
    plot sequence logo for test results
}
\description{
    plot sequence logo for test results
}
\usage{
dagLogo(testDAUresults, type=c("diff", "zscore"), pvalueCutoff=0.05, namehash=NULL,
                    font="Helvetica-Bold", textgp=gpar(), legend=FALSE, 
                    labelRelativeToAnchor=FALSE)
}
\arguments{
  \item{testDAUresults}{output of \code{\link{testDAU}}, should be an object of testDAUresults}
  \item{type}{"diff" or "zscore"}
  \item{pvalueCutoff}{pvalue cutoff for logo plot}
  \item{namehash}{the hash table to convert rownames of test results to a single letter to be plotted in the logo}
  \item{font}{font for logo symbol}
  \item{textgp}{text parameter}
  \item{legend}{plot legend or not, default false.}
  \item{labelRelativeToAnchor}{plot label relative to anchor or not, default false.}
}
\value{
    none
}
\author{
    Jianhong Ou
}
\seealso{
\code{\link{nameHash}} 
}
\examples{
    data("seq.example")
    data("proteome.example")
    bg <- buildBackgroundModel(seq.example, proteome=proteome.example, permutationSize=10)
    t <- testDAU(seq.example, bg)
    dagLogo(t)
}
\keyword{ figure }