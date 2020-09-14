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
#' Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. \emph{Biological Journal of the Linnean Society}, \bold{118}, 131-151.
#'
#' @examples
#'
#' # Get morphological distances for Michaux (1989) data set:
#' distances <- calculate_morphological_distances(cladistic_matrix = michaux_1989)
#'
#' # Show distances:
#' distances
#' @exportPattern "^[[:alpha:]]+"
#' @import ape
#' @import phytools
#' @import strap


#' @importFrom clipr write_clip
#' @importFrom geoscale geoscalePlot
#' @importFrom graphics layout legend lines par plot points polygon text
#' @importFrom grDevices adjustcolor chull hcl.colors rgb
#' @importFrom methods hasArg
# @importFrom rgl plot3d lines3d points3d text3d view3d # Still breaks Claddis so moved to suggests for now
#' @importFrom stats as.dist dist dpois pchisq runif var
#' @importFrom utils combn
NULL

#' Character-taxon matrix from Day et al. 2016
#'
#' The character-taxon matrix from Day et al. (2016).
#'
#' @name day_2016
#' @docType data
#' @format A character-taxon matrix in the format imported by \link{read_nexus_matrix}.
#' @references Day, M. O., Rubidge, B. S. and Abdala, F., 2016. A new mid-Permian burnetiamorph therapsid from the Main Karoo Basin of South Africa and a phylogenetic review of Burnetiamorpha. \emph{Acta Palaeontologica Polonica}, \bold{61}, 701-719.
#' @keywords datasets
NULL

#' Character-taxon matrix from Gauthier 1986
#'
#' The character-taxon matrix from Gauthier (1986).
#'
#' @name gauthier_1986
#' @docType data
#' @format A character-taxon matrix in the format imported by \link{read_nexus_matrix}.
#' @references Gauthier, J. A., 1986. Saurischian monophyly and the origin of birds. In Padian, K. (ed.) \emph{The Origin of Birds and the Evolution of Flight}. Towne and Bacon, San Francisco, CA, United States, 1-55.
#' @keywords datasets
NULL

#' Character-taxon matrix from Michaux 1989
#'
#' The character-taxon matrix from Michaux (1989).
#'
#' @name michaux_1989
#' @docType data
#' @format A character-taxon matrix in the format imported by
#' \link{read_nexus_matrix}.
#' @references Michaux, B., 1989. Cladograms can reconstruct phylogenies: an example from the fossil record. \emph{Alcheringa}, \bold{13}, 21-36.
#' @keywords datasets
NULL
