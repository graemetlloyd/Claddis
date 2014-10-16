\name{FindAncestor}

\alias{FindAncestor}

\title{
Find ancestor
}

\description{
Finds last common ancestor (node) of a set of descendant tips.
}

\usage{
FindAncestor(descs, tree)
}

\arguments{
\item{descs}{A vector of mode character representing the tip names for which an ancestor is sought.}
\item{tree}{The tree as a phylo object.}
}

\details{
An internal function for \link{TrimMorphDistMatrix}.
}

\value{
\item{anc.node}{The ancestral node number.}
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Nothing yet
}

\keyword{ancestor}