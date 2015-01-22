\name{AncStateEstMatrix}

\alias{AncStateEstMatrix}

\title{
Ancestral State Estimation
}

\description{
Given a tree and a cladistic matrix uses likelihood to estimate the ancestral states for every character.
}

\usage{
AncStateEstMatrix(morph.matrix, tree, estimate.allchars=FALSE, estimate.tips=FALSE)
}

\arguments{
\item{morph.matrix}{A character-taxon matrix in the format imported by \link{ReadMorphNexus}.}
\item{tree}{A tree (phylo object) with branch lengths that represents the relationships of the taxa in \code{morph.matrix}.}
\item{estimate.allchars}{An optional that allows the user to make estimates for all ancestral values. The default will only make estimates for nodes that link coded terminals.}
\item{estimate.tips}{An optional that allows the user to make estimates for tip values. The default only makes estimates for internal nodes.}
}

\details{
Uses the \link{rerootingMethod} (REF) in the \link{phytools} package.
}

\value{
\item{anc.lik.matrix}{A matrix of nodes (hypothetical ancestors; rows) against characters (columns) listing the reconstructed ancestral states.}
}

\references{
REROOTING PAPER
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Nothing yet
}
