% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.timeBins.r
\name{print.timeBins}
\alias{print.timeBins}
\title{Compact display of time bins}
\usage{
\method{print}{timeBins}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"timeBins"}.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
Nothing is directly returned, instead a text summary describing a \code{"timeBins"} object is printed to the console.
}
\description{
Displays a compact summary of a timeBins object.
}
\details{
Displays some basic summary information on a time bins object, including number of bins and their names and timespans.
}
\examples{

# Create a time bins object:
time_bins <- matrix(
  data = c(99.6, 93.5, 93.5, 89.3, 89.3, 85.8, 85.8, 83.5, 83.5, 70.6, 70.6, 65.5),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(
    c("Cenomanian", "Turonian", "Coniacian", "Santonian", "Campanian", "Maastrichtian"),
    c("fad", "lad")
  )
)

# Set class as timeBins:
class(time_bins) <- "timeBins"

# Show print.timeBins version of each included data sets:
print.timeBins(x = time_bins)
}
\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}
