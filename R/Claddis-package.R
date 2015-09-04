#' Measuring Morphological Diversity and Evolutionary Tempo
#' 
#' Measures morphological diversity from discrete character data and estimates evolutionary tempo on phylogenetic trees.
#' 
#' @name Claddis-package
#'
#' @aliases Claddis
#'
#' @docType package
#'
#' @author Graeme T. Lloyd <graemetlloyd@@gmail.com>
#'
#' @references Brusatte, S. L., Benton, M. J., Ruta, M. and Lloyd, G. T.,
#' 2008a. Superiority, competition, and opportunism in the evolutionary
#' radiation of dinosaurs. Science, 321, 1485-1488.
#' 
#' Brusatte, S. L., Benton, M. J., Ruta, M. and Lloyd, G. T., 2008b. The first
#' 50 million years of dinosaur evolution: macroevolutionary pattern and
#' morphological disparity. Biology Letters, 4, 733-736.
#' 
#' Brusatte, S. L., Lloyd, G. T., Wang, S. C. and Norell, M. A., 2014. Gradual
#' assembly of avian body plan culminated in rapid rates of evolution across
#' dinosaur-bird transition. Current Biology, 24, 2386-2392.
#'
#' Close, R. A., Friedman, M., Lloyd, G. T. and Benson, R. B. J., 2015.
#' Evidence for a mid-Jurassic adaptive radiation in mammals. Current Biology,
#' 25, 2137-2142.
#'
#' Lloyd, G. T., Wang, S. C. and Brusatte, S. L., 2012. Identifying
#' heterogeneity in rates of morphological evolution: discrete character change
#' in the evolution of lungfish (Sarcopterygii; Dipnoi). Evolution, 66,
#' 330-348.
#'
#' @keywords disparity,distance,morphology,phylogeny
#'
#' @examples
#' 
#' # Get morphological distances for Michaux (1989) data set:
#' distances <- MorphDistMatrix(Michaux1989)
#' 
#' # Show distances:
#' distances
#'
#' @exportPattern "^[[:alpha:]]+"
#' @import ape
#' @import gdata
#' @import phytools
#' @import strap
NULL

#' @importFrom stats dist
NULL

#' @importFrom stats dpois
NULL

#' @importFrom stats pchisq
NULL

#' @importFrom stats runif
NULL

#' @importFrom utils combn
NULL


#' Character-taxon matrix from Gauthier 1986
#' 
#' The character-taxon matrix from Gauthier (1986).
#' 
#' 
#' @name Gauthier1986
#' @docType data
#' @format A character-taxon matrix in the format imported by \link{ReadMorphNexus}.
#' @references Gauthier, J. A., 1986. Saurischian monophyly and the origin of
#' birds. In Padian, K. (ed.) The Origin of Birds and the Evolution of Flight.
#' Towne and Bacon, San Francisco, CA, United States, 1-55.
#' @keywords datasets
NULL





#' Character-taxon matrix from Michaux 1989
#' 
#' The character-taxon matrix from Michaux (1989).
#' 
#' 
#' @name Michaux1989
#' @docType data
#' @format A character-taxon matrix in the format imported by
#' \link{ReadMorphNexus}.
#' @references Michaux, B., 1989. Cladograms can reconstruct phylogenies: an
#' example from the fossil record. Alcheringa, 13, 21-36.
#' @keywords datasets
NULL



