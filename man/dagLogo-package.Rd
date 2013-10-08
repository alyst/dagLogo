\name{dagLogo-package}
\alias{dagLogo-package}
\docType{package}
\title{
Visualize significant conserved amino acid sequence pattern in groups 
based on probability theory
}
\description{
We implement iceLogo by R to visualize significant conserved amino acid sequence pattern
based on probability theory. Compare to iceLogo, dagLogo can also visualize significant 
sequence patterns by clustering the peptides by groups such as charge, chemistry, hydrophobicity
and etc.
}
\details{
\tabular{ll}{
Package: \tab dagLogo\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-09-31\cr
License: \tab GPL (>= 2)\cr
}
DAG: Differential Amino acid Group

There are several differences between dagLogo from iceLogo:

1. The sequence patterns can be grouped by charge, chemistry, hydrophobicity and etc.

2. dagLogo accepts different length of aligned amino acid sequences.

3. Except Random, regional (called restricted in dagLogo) and terminal 
(called anchored) background model, the background sequence could be set to 
other regions of the genes in inputs and complementary set of the proteome.
}
\author{
Jianhong Ou, Julie Lihua Zhu

Maintainer: Jianhong Ou <jianhong.ou@umassmed.edu>
}
\keyword{ package }
\examples{
    data("seq.example")
    data("proteome.example")
    bg <- buildBackgroundModel(seq.example, proteome=proteome.example, permutationSize=10L)
    t <- testDAU(seq.example, bg)
    dagLogo(t)
}
