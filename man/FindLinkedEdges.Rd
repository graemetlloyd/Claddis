\name{FindLinkedEdges}

\alias{FindLinkedEdges}

\title{
Find linked edges for a tree
}

\description{
Given a tree finds edges that are linked to each other.
}

\usage{
FindLinkedEdges(tree)
}

\arguments{
\item{tree}{A tree (phylo object).}
}

\details{
Finds all edges that link (share a node) with each edge of a tree.

This is intended as an internal function, but may be of use to someone else.
}

\value{
Returns a matrix where links are scored 1 and everything else 0. The diagonal is left as zero.
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Create a simple four-taxon tree:
tree <- read.tree(text="(A,(B,(C,D)));")

# Find linked (1) edges matrix for tree:
FindLinkedEdges(tree)
}
