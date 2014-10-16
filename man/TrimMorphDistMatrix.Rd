\name{TrimMorphDistMatrix}

\alias{TrimMorphDistMatrix}

\title{
Trim morphological distance matrix
}

\description{
Trims a morphological distance matrix by removing nodes that cause empty cells.
}

\usage{
TrimMorphDistMatrix(dist.matrix, tree=NULL)
}

\arguments{
\item{dist.matrix}{A distance matrix in the format created by \link{MorphDistMatrix}.}
\item{tree}{If the distance matrix includes ancestors this should be the tree (phylo object) used to reconstruct them.}
}

\details{
Trims a morphological distance matrix by removing nodes that cause empty cells.

Some distances are not calculable from cladistic matrices if there are taxa that have no coded characters in common. This algorithm iteratively removes the taxa responsible for the most empty cells until the matrix is complete.

If the matrix includes reconstructed ancestors the user should also provide the tree used (as the \code{tree} argument). The function will then also remove the tips from the tree and where reconstructed ancestors also cause empty cells will prune the minimum number of taxa required to also remove these. The function will then renumber the nodes in the distance matrix so they match the pruned tree.
}

\value{
\item{dist.matrix}{A complete distance matrix with all cells filled. If there were no empty cells will return.}
\item{tree}{A tree (if supplied) with the removed taxa (see below) pruned. If no taxa are dropped will return the same tree as inputted. If no tree is supplied this is set to NULL.}
\item{removed.taxa}{A vector of mode character listing the taxa removed. If none are removed this will be set to NULL.}
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Nothing yet
}

\keyword{distance}