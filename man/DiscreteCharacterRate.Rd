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
\item{tree}{A tree (phylo object) with branch lengths that represents the relationships of the taxa in \code{morph.matrix}.}
\item{clad.matrix}{A character-taxon matrix in the format imported by \link{ReadMorphNexus}.}
\item{time.bins}{An optional that allows the user to make estimates for all ancestral values. The default will only make estimates for nodes that link coded terminals.}
\item{alpha}{The alpha value to be used for the significance tests. The default is 0.01.}
}

\details{
\link{AncStateEstMatrix}
}

\value{
\item{node.results}{Text.}
\item{branch.results}{Text.}
\item{per.bin.rates}{Text.}
}

\references{
Lloyd et al. 2012
Brusatte et al. in press
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com} and Steve C. Wang \email{scwang@swarthmore.edu}
}

\examples{
# Nothing yet
}

\keyword{evolutionary rates}