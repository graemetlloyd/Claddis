# devtools::build_win(version = "R-devel")

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
#' @references
#'
#' Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. Biological Journal of the Linnean Society, 118, 131-151.
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
#' @import phytools
#' @import strap


# there's no hope, move rgl to suggests!
# @import rgl
# @importFrom rgl plot3d lines3d points3d text3d view3d 

#' @importFrom clipr write_clip
#' @importFrom gdata trim
#' @importFrom graphics layout lines par plot points polygon text 
#' @importFrom grDevices adjustcolor chull 
#' @importFrom stats as.dist dist dpois pchisq runif var 
#' @importFrom utils combn 
NULL

#' Character-taxon matrix from Day et al. 2016
#'
#' The character-taxon matrix from Day et al. (2016).
#'
#'
#' @name Day2016
#' @docType data
#' @format A character-taxon matrix in the format imported by \link{ReadMorphNexus}.
#' @references Day, M. O., Rubidge, B. S. and Abdala, F., 2016. A new mid-Permian
#' burnetiamorph therapsid from the Main Karoo Basin of South Africa and a
#' phylogenetic review of Burnetiamorpha. Acta Palaeontologica Polonica, 61, 701-719.
#' @keywords datasets
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



