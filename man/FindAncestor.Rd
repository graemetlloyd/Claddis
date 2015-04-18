\name{FindAncestor}

\alias{FindAncestor}

\title{
Find ancestor
}

\description{
Finds last common ancestor (node) of a set of two or more descendant tips.
}

\usage{
FindAncestor(descs, tree)
}

\arguments{
\item{descs}{A vector of mode character representing the tip names for which an ancestor is sought.}
\item{tree}{The tree as a phylo object.}
}

\details{
Intended for use as an internal function for \link{TrimMorphDistMatrix}, but potentially of more general use.
}

\value{
\item{anc.node}{The ancestral node number.}
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Create a simple four-taxon tree:
tree <- read.tree(text="(A,(B,(C,D)));")

# Plot the tree:
plot(tree)

# Add nodelabels and show that the most recent common
# ancestor of B, C, and D is node 6:
nodelabels()

# Use FindAncestor to show that the most recent common
# ancestor of B, C, and D is node 6:
FindAncestor(c("B", "C", "D"), tree)
}

\keyword{ ancestor }
