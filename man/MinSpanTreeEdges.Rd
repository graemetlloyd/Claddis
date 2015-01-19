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
# Nothing yet
}
