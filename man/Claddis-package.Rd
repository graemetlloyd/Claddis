\name{Claddis-package}

\alias{Claddis-package}

\alias{Claddis}

\docType{package}

\title{
Cladistic disparity
}

\description{
A package to read cladistic morphological datasets into R and 
perform morphological disparity or rates of character 
evolution analyses upon them.
}

\details{
\tabular{ll}{
Package: \tab Claddis\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2013-05-22\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
}

\author{
Graeme T. Lloyd <graemetlloyd@gmail.com>
}

\references{
Brusatte, S. L., Benton, M. J., Ruta, M. and Lloyd, G. T., 2008a. Superiority, competition, and opportunism in the evolutionary radiation of dinosaurs. Science, 321, 1485-1488.

Brusatte, S. L., Benton, M. J., Ruta, M. and Lloyd, G. T., 2008b. The first 50 million years of dinosaur evolution: macroevolutionary pattern and morphological disparity. Biology Letters, 4, 733-736.

Brusatte, S. L., Lloyd, G. T., Wang, S. C. and Norell, M. A., 2014. Gradual assembly of avian body plan culminated in rapid rates of evolution across dinosaur-bird transition. Current Biology, 24, 2386-2392.

Lloyd, G. T., Wang, S. C. and Brusatte, S. L., 2012. Identifying heterogeneity in rates of morphological evolution: discrete character change in the evolution of lungfish (Sarcopterygii; Dipnoi). Evolution, 66, 330-348.
}

\keyword{ disparity,distance,morphology,phylogeny }

\examples{
# Get morphological distances for Michaux (1989) data set:
distances <- MorphDistMatrix(Michaux1989)

# Show distances:
distances
}
