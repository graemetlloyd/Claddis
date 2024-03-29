% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Claddis-package.R
\docType{package}
\name{Claddis-package}
\alias{Claddis-package}
\alias{Claddis}
\title{Measuring Morphological Diversity and Evolutionary Tempo
\if{html}{\figure{logo.png}{options: style='float: right' alt='logo' width='120'}}
Measures morphological diversity from discrete character data and estimates evolutionary tempo on phylogenetic trees.}
\examples{

# Get morphological distances for Michaux (1989) data set:
distances <- calculate_morphological_distances(cladistic_matrix = michaux_1989)

# Show distances:
distances
}
\references{
Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. \emph{Biological Journal of the Linnean Society}, \bold{118}, 131-151.
}
\author{
Graeme T. Lloyd <graemetlloyd@gmail.com>
}
