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
Uses \link{AncStateEstMatrix}.
}

\value{
\item{node.results}{A table displaying the results of the per-clade rate tests.}
\item{branch.results}{A table displaying the results of the per-branch rate tests.}
\item{per.bin.rates}{In preparation. DO NOT USE.}
}

\references{
Brusatte et al. 2014
Lloyd et al. 2012
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com} and Steve C. Wang \email{scwang@swarthmore.edu}
}

\examples{
# Nothing yet
}

\keyword{evolutionary rates}