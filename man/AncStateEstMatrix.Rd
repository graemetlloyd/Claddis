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
Uses the \link{rerootingMethod} (Tang et al. 1995) as implemented in the \link{phytools} package to make ancestral state estimates. Here these are collapsed to the most likely state, or if two or more states are most likely, a polymorphism of the most likely states. This is the method used by Brusatte et al. (2014).
}

\value{
\item{anc.lik.matrix}{A matrix of nodes (hypothetical ancestors; rows) against characters (columns) listing the reconstructed ancestral states.}
}

\references{
Brusatte, S. L., Lloyd, G. T., Wang, S. C. and Norell, M. A., 2014. Gradual assembly of avian body plan culminated in rapid rates of evolution across dinosaur-bird transition. Current Biology, 24, 2386-2392.

Yang, Z., Kumar, S. and Nei, M., 1995. A new method of inference of ancestral nucleotide and amino acid sequences. Genetics, 141, 1641-1650.
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Set random seed:
set.seed(17)

# Generate a random tree for the Michaux data set:
tree <- rtree(nrow(Michaux1989$matrix))

# Update taxon names to match those in the data matrix:
tree$tip.label <- rownames(Michaux1989$matrix)

# Set root time by making youngest taxon extant:
tree$root.time <- max(diag(vcv(tree)))

# Estimate ancestral states:
AncStateEstMatrix(Michaux1989, tree)
}
