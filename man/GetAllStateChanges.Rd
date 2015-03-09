\name{GetAllStateChanges}

\alias{GetAllStateChanges}

\title{
Finds all state changes on a tree using stochastic character mapping
}

\description{
Takes a cladistic matrix and time-scaled tree and makes points estimates for every character change using stochastic character mapping.
}

\usage{
GetAllStateChanges(clad.matrix, tree, time.bins, Nsim=10)
}

\arguments{
\item{clad.matrix}{A character-taxon matrix in the format imported by \link{ReadMorphNexus}.}
\item{tree}{A time-scaled tree (phylo object) that represents the relationships of the taxa in \code{clad.matrix}.}
\item{time.bins}{A vector of ages representing the boundaries of a series of time bins.}
\item{Nsim}{The number of simulations to perform (passed to \code{make.simmap}.}
}

\details{
A wrapper function for \link{make.simmap} in the \link{phytools} package.

More to come.
}

\value{
\item{all.state.changes}{A matrix of rows for each change with columns corresponding to: the character concerned, the simulation number, the edge number, the time, and the start and end states.}
\item{character.times}{A vector of the sampled tree-length (in Ma) for each character.}
\item{edge.length.per.bin}{A matrix of time bins (columns) and characters (rows) indicating the sampled tree-length (in Ma).}
\item{terminal.edge.length.per.bin}{As above, but for terminal edges only.}
\item{internal.edge.length.per.bin}{As above, but for internal edges only.}
}

\references{
STOCHASTIC CHARACTER MAPPING
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Nothing yet
}
