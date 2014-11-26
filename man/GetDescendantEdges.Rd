\name{GetDescendantEdges}

\alias{GetDescendantEdges}

\title{
Gets descendant edges of an internal node
}

\description{
Gets all descendant edges of a phylo object for an internal node.
}

\usage{
GetDescendantEdges(n, tree)
}

\arguments{
\item{n}{An integer corresponding to the internal node for which the descendant edges are sought.}
\item{tree}{A tree as a phylo object.}
}

\details{
Returns a vector of integers corresponding to row numbers in \code{$edge} or cells in \code{$edge.length} of the descendant edges of the internal node supplied.
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Nothing yet
}