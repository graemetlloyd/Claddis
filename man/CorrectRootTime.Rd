\name{CorrectRootTime}

\alias{CorrectRootTime}

\title{
Corrects root.time after taxa have been pruned from a tree
}

\description{
Corrects root.time after taxa have been pruned from a tree using drop.tip
}

\usage{
CorrectRootTime(original.tree, pruned.tree)
}

\arguments{
\item{original.tree}{A tree in phylo format.}
\item{pruned.tree}{A tree in phylo format that represents a pruned version of \code{original.tree}.}
}

\details{
When removing taxa from a time-scaled tree using \link{drop.tip} in \link{ape} \code{$root.time} is left unchanged. This can cause downstream problems if not corrected and that is what this function does.

Note that \code{fixRootTime} in the \code{paleotree} package performs the same function, but is not called here to reduce the number of libraries on which \code{Claddis} is dependent.
}

\value{
Returns a tree (phylo object) with a corrected \code{$root.time}.
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Nothing yet
}
