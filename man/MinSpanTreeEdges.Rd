\name{MinSpanTreeEdges}

\alias{MinSpanTreeEdges}

\title{
Get edges of minimum spanning tree
}

\description{
Returns edges of a minimum spanning tree given a distance matrix.
}

\usage{
MinSpanTreeEdges(dist.matrix)
}

\arguments{
\item{dist.matrix}{A square matrix of distances between objects.}
}

\details{
This function is a wrapper for \link{mst} in the \link{ape} package, but returns a vector of edges rather than a square matrix of links.
}

\value{
A vector of named edges (X->Y) with their distances. The sum of this vector is the length of the minimum spanning tree.
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Create a simple square matrix of distances:
dist.matrix <- matrix(c(0,1,2,3,1,0,1,2,2,1,0,1,3,2,1,0), nrow=4,
  dimnames=list(LETTERS[1:4], LETTERS[1:4]))

# Show matrix to confirm that the off diagonal has the shortest
# distances:
dist.matrix

# Use MinSpanTreeEdges to get the edges for the minimum spanning
# tree:
MinSpanTreeEdges(dist.matrix)

# Use sum of MinSpanTreeEdges to get the length of the minimum
# spanning tree:
sum(MinSpanTreeEdges(dist.matrix))
}
