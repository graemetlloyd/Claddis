\name{ChangesInBins}

\alias{ChangesInBins}

\title{
Counts the changes in a series of time bins
}

\description{
Given a vector of dates for a series of time bins and another for the times when a character change occurred will returns the total number of changes in each bin.
}

\usage{
ChangesInBins(change.times, time.bins)
}

\arguments{
\item{change.times}{A vector of ages in millions of years at which character changes are hypothesised to have occurred.}
\item{time.bins}{A vector of ages in millions of years of time bin boundaries in old-to-young order.}
}

\details{
Calculates the total number of evolutionary changes in a series of time bins. This is intended as an internal function for rate calculations, but may be of use to someone.
}

\value{
A vector giving the number of changes for each time bin. Names indicate the maximum and minimum values for each time bin.
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\examples{
# Create a random dataset of 100 changes:
change.times <- runif(100, 0, 100)

# Create time bins:
time.bins <- seq(100, 0, length.out=11)

# Get N changes for each bin:
ChangesInBins(change.times, time.bins)
}
