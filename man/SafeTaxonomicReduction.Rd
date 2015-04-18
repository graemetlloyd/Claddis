\name{SafeTaxonomicReduction}

\alias{SafeTaxonomicReduction}

\title{
Safe Taxonomic Reduction
}

\description{
Performs Safe Taxonomic Reduction (STR) on a character-taxon matrix.
}

\usage{
SafeTaxonomicReduction(morph.matrix)
}

\arguments{
\item{morph.matrix}{A character-taxon matrix in the format imported by \link{ReadMorphNexus}.}
}

\details{
Performs Safe Taxonomic Reduction (Wilkinson 1995).

If no taxa can be safely removed will simply return the text: "No taxa can be safely removed".
}

\value{
\item{str.list}{A matrix listing the taxa that can be removed (\code{Junior}), the taxa which they are equivalent to (\code{Senior}) and the rule under which they can be safely removed (\code{Rule}).}
\item{reduced.matrix}{A character-taxon matrix excluding the taxa that can be safely removed.}
\item{removed.matrix}{A character-taxon matrix of the taxa that can be safely removed.}
}

\references{
Wilkinson, M., 1995. Coping with abundant missing entries in phylogenetic inference using parsimony. Systematic Biology, 44, 501-514.
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Performs STR on the Gauthier 1986 dataset used
# in the Wilkinson (1995) paper:
str.out <- SafeTaxonomicReduction(Gauthier1986)

# View deleted taxa:
str.out$str.list

# View reduced matrix:
str.out$reduced.matrix

# View removed matrix:
str.out$removed.matrix
}

\keyword{Safe Taxonomic Reduction}
