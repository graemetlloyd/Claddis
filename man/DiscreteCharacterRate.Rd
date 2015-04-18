\name{DiscreteCharacterRate}

\alias{DiscreteCharacterRate}

\title{
A function for examining discrete character rates
}

\description{
Given a tree and a cladistic matrix uses likelihood methods to search for significant rate excursions from an equal rate model.
}

\usage{
DiscreteCharacterRate(tree, clad.matrix, time.bins, alpha=0.01)
}

\arguments{
\item{tree}{A tree (phylo object) with branch lengths that represents the relationships of the taxa in \code{clad.matrix}.}
\item{clad.matrix}{A character-taxon matrix in the format imported by \link{ReadMorphNexus}.}
\item{time.bins}{A vector of ages (in millions of years) indicating the boundaries of a series of time bins.}
\item{alpha}{The alpha value to be used for the significance tests. The default is 0.01.}
}

\details{
Uses \link{AncStateEstMatrix} to estimate ancestral states and parsimony to count changes on branches, then applies the methods outlined in Lloyd et al. (2012), and updated in Brusatte et al. (2014) to test for significantly high or low rates on a branch or within a clade.

Please note that the time series method used here is an untested attempt to deal with the problem of correcting for completeness (something that Lloyd et al. 2012 were not able to do) and should be used with caution.
}

\value{
\item{node.results}{A table displaying the results of the per-clade rate tests.}
\item{branch.results}{A table displaying the results of the per-branch rate tests.}
\item{per.bin.rates}{Per time-bin rates (use with caution).}
}

\references{
Brusatte, S. L., Lloyd, G. T., Wang, S. C. and Norell, M. A., 2014. Gradual assembly of avian body plan culminated in rapid rates of evolution across dinosaur-bird transition. Current Biology, 24, 2386-2392.

Lloyd, G. T., Wang, S. C. and Brusatte, S. L., 2012. Identifying heterogeneity in rates of morphological evolution: discrete character change in the evolution of lungfish (Sarcopterygii; Dipnoi). Evolution, 66, 330-348.
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com} and Steve C. Wang \email{scwang@swarthmore.edu}
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

# Get discrete character rates:
DiscreteCharacterRate(tree, Michaux1989,
  seq(tree$root.time, 0, length.out=3), alpha=0.01)
}

\keyword{ evolution,rates }
