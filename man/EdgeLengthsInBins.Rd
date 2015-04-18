\name{EdgeLengthsInBins}

\alias{EdgeLengthsInBins}

\title{
Edge-lengths present in time-bins
}

\description{
Given a time-scaled tree and set of time bin boundaries will sum the edge-lengths present in each bin.
}

\usage{
EdgeLengthsInBins(tree, time.bins, pruned.tree=NULL)
}

\arguments{
\item{tree}{A time-scaled tree in phylo format with a \code{$root.time} value.}
\item{time.bins}{A vector of ages in millions of years of time bin boundaries in old-to-young order.}
\item{pruned.tree}{A time-scaled tree in phylo format with a \code{$root.time} value that is a subset of \code{tree}.}
}

\details{
Calculates the total edge length present in each of a series of time bins. This is intended as an internal function for rate calculations, but may be of use to someone.

The option of using a \code{pruned.tree} allows the user to correctly classify internal and terminal branches in a subtree of the larger tree. So for example, if taxa A and B are sisters then after pruning B the subtree branch leading to A is composed of an internal and a terminal branch on the complete tree. 
}

\value{
\item{edge.length.in.bin}{A vector giving the summed values in millions of years for each time bin. Names indicate the maximum and minimum values for each time bin.}
\item{terminal.edge.length.in.bin}{As above, but counting terminal edges only.}
\item{internal.edge.length.in.bin}{As above, but counting internal edges only.}
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Create a random 10-taxon tree:
tree <- rtree(10)

# Add root age:
tree$root.time <- 100

# Create time bins:
time.bins <- seq(100, 0, length.out=11)

# Get edge lengths for each bin:
EdgeLengthsInBins(tree, time.bins)
}
