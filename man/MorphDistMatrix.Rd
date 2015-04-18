\name{MorphDistMatrix}

\alias{MorphDistMatrix}

\title{
Get distance matrices from a cladistic matrix
}

\description{
Takes a cladistic morphological dataset and converts it into a set of pairwise distances.
}

\usage{
MorphDistMatrix(morph.matrix, transform.proportional.distances="arcsine_sqrt")
}

\arguments{
\item{morph.matrix}{A character-taxon matrix in the format imported by \link{ReadMorphNexus}.}
\item{transform.proportional.distances}{Whether to transform the proportional distances (Gower and max). Options are \code{none}, \code{sqrt}, or \code{arcsine_sqrt} (the default).}
}

\details{
This function is an important preliminary step in performing multivariate ordination(s) upon a cladistic dataset of discrete characters and deals with three peculiarities of such data: 1) the prevalence of missing values, 2) the use of unordered multistate characters, and 3) the presence of polymorphisms.

Missing data is dealt with by providing three rescaled distances as well as the uncorrected raw distances. These are: 1) the Generalised Euclidean Distance (GED) of Wills (2001) which replaces all incalculable taxon-taxon character distances with a weighted fractional mean distance, 2) Gower (1971) dissimilarity, which rescales each pairwise distance by dividing by the number of comparable characters upon which that distance is based, and 3) a set of distances rescaled against the maximum possible observable distance based on the set of comparable characters upon which that distance is based (using the difference between \code{max.vals} and \code{min.vals} supplied by the input matrix). In practice the user is not recommended to use raw distances unless the input matrix is complete.

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

Wills, M. A., 2001a. Morphological disparity: a primer. In: Adrain, J. M., Edgecombe, G. D. and Lieberman, B. S. (eds.), Fossils, Phylogeny, and Form: An Analytical Approach. Kluwer Academic/Plenum Publishers, New York, p55-144.
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Get morphological distances for Michaux (1989)
# data set:
distances <- MorphDistMatrix(Michaux1989)

# Show raw Euclidean distances:
distances$raw.dist.matrix

# Show Generailsed Euclidean Distances:
distances$GED.dist.matrix

# Show Gower rescaled distances:
distances$gower.dist.matrix

# Show maximum observable rescaled distances:
distances$max.dist.matrix

# Show number of characters that can be scored for
# each pairwise comparison:
distances$comp.char.matrix
}

\keyword{ distance }
