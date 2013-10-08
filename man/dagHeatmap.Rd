\name{dagHeatmap}
\alias{dagHeatmap}
\title{
    plot heatmap for test results
}
\description{
    plot heatmap for test results
}
\usage{
dagHeatmap(testDAUresults, type=c("diff", "zscore"), \dots)
}
\arguments{
  \item{testDAUresults}{output of \code{\link{testDAU}}, should be an object of testDAUresults}
  \item{type}{"diff" or "zscore"}
  \item{\dots}{parameter could be passed to pheatmap}
}
\value{
    none
}
\author{
    Jianhong Ou
}
\examples{
    data("seq.example")
    data("proteome.example")
    bg <- buildBackgroundModel(seq.example, proteome=proteome.example, permutationSize=10)
    t <- testDAU(seq.example, bg)
    dagHeatmap(t)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ figure }