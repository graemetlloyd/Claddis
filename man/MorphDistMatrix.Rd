\name{MorphDistMatrix}

\alias{MorphDistMatrix}

\title{
Get distance matrices from a character-taxon matrix
}

\description{
Takes a cladistic morphological dataset of discrete characters and converts it into a set of pairwise distances.
}

\usage{
MorphDistMatrix(morph.matrix)
}

\arguments{
\item{morph.matrix}{A character-taxon matrix in the format imported by \link{ReadMorphNexus}.}
}

\details{
This function is an important preliminary step in performing multivariate ordination(s) upon a cladistic dataset of discrete characters and deals with three peculiarities of such data: 1) the prevalence of missing values, 2) the use of unordered multistate characters, and 3) the presence of polymorphisms.

Missing data is dealt with by providing two rescaled distances as well as the uncorrected raw distances. These are: 1) Gower (1971) dissimilarity, which rescales each pairwise distance by dividing by the number of comparable characters upon which that distance is based, and 2) a set of distances rescaled against the maximum possible observed distance based on the set of comparable characters upon which that distance is based (using the difference between \code{max.vals} and \code{min.vals} supplied by the input matrix). In practice the user is not recommended to use raw distances unless the input matrix is highly complete.

The method automatically treats distances between unordered characters as zero (if the two states are identical) or one (if the two states are different). So, for example, the distances between 0 and 3 and between 2 and 3 for an unordered character are both 1.

Finally, polymorphisms are dealt with by using the minimum possible distance by considering all possible values implied by the polymorphism. So, for example, the distance between (01) and 3 will be considered to be 1 for an unordered character and 2 (the minimum distance, between 1 and 3) for an ordered character.

All metrics are rescaled according to character weightings, as supplied by the \code{weights} vector of the input matrix.

It is important to note that in practice not all pairwise distances can be calculated due to missing data. In these cases incomplete distance matrices will be returned, with incalculable values scored as \code{NA}. In such cases the user is advised to apply the \link{TrimMorphDistMatrix} function before ordination.
}

\value{
\item{raw.dist.matrix}{A distance matrix indicating the raw distances (based on the method set in \code{dist.method}) of each pairwise comparison.}
\item{GED.dist.matrix}{A Generalised Euclidean Distance (GED) matrix (Wills 2001).}
\item{gower.dist.matrix}{A distance matrix where raw distances have been rescaled using Gower (1971) dissimilarity (then rescaled from 0 to 1 and the arcsine of the square root taken).}
\item{max.dist.matrix}{A distance matrix where raw distances have been rescaled against the maximum possible observable distance (then the arcsine of the square root taken).}
\item{comp.char.matrix}{A matrix showing the number of coded characters that were used to make each pairwise comparison.}
}

\references{
Gower, J. C., 1971. A general coefficient of similarity and some of its properties. Biometrics, 27, 857-871.
Wills, M. A., 2001.
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Nothing yet
}

\keyword{Distance}