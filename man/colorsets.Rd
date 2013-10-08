\name{colorsets}
\alias{colorsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
retrieve color setting for logo
}
\description{
retrieve prepared color setting for logo
}
\usage{
colorsets(colorScheme=c("null", "classic", "charge", "chemistry", "hydrophobicity"))
}
\arguments{
  \item{colorScheme}{could be 'null', 'charge', 'chemistry', 'classic' or 'hydrophobicity'}
}
\value{
    A character vector of color scheme
}
\author{
Jianhong Ou
}
\examples{
col <- colorsets("hydrophobicity")
}
\keyword{figure}