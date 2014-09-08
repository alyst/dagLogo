\name{nameHash}
\alias{nameHash}
\title{
convert group name to a single character
}
\description{
convert group name to a single character to shown in a logo
}
\usage{
nameHash(nameScheme=c("classic", "charge", "chemistry", "hydrophobicity"))
}
\arguments{
  \item{nameScheme}{could be "classic", "charge", "chemistry", "hydrophobicity"}
}
\value{
    A character vector of name scheme
}
\author{
    Jianhong Ou
}
\examples{
    nameHash("charge")
}
\keyword{ figure }