\name{GetNodeAges}

\alias{GetNodeAges}

\title{
Returns node ages for a time-scaled tree
}

\description{
Given a tree with branch-lengths scaled to time and a value for \code{$root.time} will return a vector of node ages.
}

\usage{
GetNodeAges(tree)
}

\arguments{
\item{tree}{A tree (phylo object) with branch lengths representing time and a value for \code{$root.time}.}
}

\details{
Returns a vector of node ages (terminal and internal) labelled by their node number.
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Create simple four-taxon tree with edge lengths all
# set to 1 Ma:
tree <- read.tree(text="(A:1,(B:1,(C:1,D:1):1):1);")

# Set root.time as 10 Ma:
tree$root.time <- 10

# Get node ages:
GetNodeAges(tree)
}
