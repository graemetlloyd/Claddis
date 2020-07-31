#' Principal Coordinates on a Cladistic Matrix
#'
#' @description
#'
#' Performs Principal Coordinates Analysis (PCoA) on a cladistic matrix.
#'
#' @param cladistic.matrix A vector of mode character representing the tip names for which an ancestor is sought.
#' @param distance.metric The distance method to use (one of "RED", "GED", "GC", or "MORD" - the default). See \link{calculate_morphological_distances} for more details.
#' @param ged.type The type of GED use. Must be one of \code{"Legacy"}, \code{"Hybrid"}, or \code{"Wills"} (the default). See details for an explanation.
#' @param distance.transformation The transformation to apply to distances. See \link{calculate_morphological_distances} for details.
#' @param distance.polymorphism.behaviour The distance behaviour for dealing with polymorphisms. Must be one of \code{"mean.difference"}, \code{"min.difference"} (the default), or \code{"random"}. See \link{calculate_morphological_distances} for details.
#' @param distance.uncertainty.behaviour The distance behaviour for dealing with uncertainties. Must be one of \code{"mean.difference"}, \code{"min.difference"} (the default), or \code{"random"}. See \link{calculate_morphological_distances} for details.
#' @param distance.inapplicable.behaviour The behaviour for dealing with inapplicables. Must be one of \code{"missing"} (default), or \code{"HSJ"}. See \link{calculate_morphological_distances} for details.
#' @param character.dependencies Only relevant if using \code{inapplicable.behaviour = "HSJ"}. Must be a two-column matrix with colnames "DependentCharacter" and "IndependentCharacter" that specifies character hierarchies. See \link{calculate_morphological_distances} for details.
#' @param alpha The alpha value (sensu Hopkins and St John 2018). Only relevant if using \code{inapplicable.behaviour = "HSJ"}. See \link{calculate_morphological_distances} for details.
#' @param correction The negative eigenvalue correction to use (one of "lingoes", "none", or "cailliez" - the default). See \link{pcoa} for more details.
#' @param time_tree If a phylmorphospace is desired then a tree with root age and branch-lengths must be included.
#' @param estimate.all.nodes If including a tree whether you want to estinate ancestral states for all characters (default is FALSE). See \link{estimate_ancestral_states} for more details.
#' @param estimate.tip.values If including a tree whether you want to estinate missing or polymorphic tip states (default is FALSE). See \link{estimate_ancestral_states} for more details.
#' @param inapplicables.as.missing See \link{estimate_ancestral_states}.
#' @param ancestral.polymorphism.behaviour Behaviour for dealing with polymorphisms when producing ancestral state estimates - see \link{estimate_ancestral_states}.
#' @param ancestral.uncertainty.behaviour Behaviour for dealing with uncertainties when producing ancestral state estimates - see \link{estimate_ancestral_states}.
#' @param threshold threshold for ancestral state estimation of discrete characters - see \link{estimate_ancestral_states} for details.
#' @param allow.all.missing Logical to allow all missing character values - see \link{estimate_ancestral_states} for details.
#'
#' @details
#'
#' Takes a cladistic matrix in the format imported by \link{read_nexus_matrix} and performs Principal Coordinates (Gower 1966) analysis on it.
#'
#' This function is effectively a wrapper for \link{pcoa} from the \link{ape} package and the user is referred there for some of the options (e.g., using the Caillez 1983 approach to avoiding negative eigenvalues).
#'
#' If providing a tree and inferring ancestral states then options to also infer missing or uncertain tips and whether to infer values for all characters at all internal nodes are provided (via \link{estimate_ancestral_states}).
#'
#' Other options within the function concern the distance metric to use and the transformation to be used if selecting a propotional distance (see \link{calculate_morphological_distances}).
#'
#' IMPORTANT: The function can remove taxa (or if including a tree, nodes as well) if they lead to an incomplete distance matrix (see \link{trim_matrix}).
#'
#' @return \item{time_tree}{The tree (if supplied). Note this may be pruned from the input tree by \link{trim_matrix}.}
#' @return \item{DistMatrix}{The distance matrix. Note this may be pruned by \link{trim_matrix} and thus not include all taxa.}
#' @return \item{RemovedTaxa}{A vector of taxa (or nodes) removed by \link{trim_matrix}. Returns NULL if none are removed.}
#' @return \item{note}{See \link{pcoa}.}
#' @return \item{values}{See \link{pcoa}.}
#' @return \item{vectors}{See \link{pcoa}.}
#' @return \item{trace}{See \link{pcoa}.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Cailliez, F., 1983. The analytical solution of the additive constant problem. \emph{Psychometrika}, \bold{48}, 305-308.
#'
#' Gower, J. C., 1966. Some distance properties of latent root and vector methods used in multivariate analysis. \emph{Biometrika}, \bold{53}, 325-338.
#'
#' @examples
#'
#' # Run on Michaux (1989) data set with defaults:
#' x <- ordinate_cladistic_matrix(Michaux1989)
#'
#' # Show output:
#' x
#'
#' # Generate a (made up) tree:
#' time_tree <- rtree(length(rownames(Michaux1989$matrix_1$matrix)))
#'
#' # Add taxon names to it:
#' time_tree$tip.label <- rownames(Michaux1989$matrix_1$matrix)
#'
#' # Set root time by making youngest taxon extant:
#' time_tree$root.time <- max(diag(vcv(time_tree)))
#'
#' # Run with tree:
#' y <- ordinate_cladistic_matrix(Michaux1989, time_tree = time_tree)
#'
#' # Show new output:
#' y
#'
#' @export ordinate_cladistic_matrix
ordinate_cladistic_matrix <- function(cladistic.matrix, distance.metric = "MORD", ged.type = "Wills", distance.transformation = "arcsine_sqrt", distance.polymorphism.behaviour = "min.difference", distance.uncertainty.behaviour = "min.difference", distance.inapplicable.behaviour = "missing", character.dependencies = NULL, alpha = 0.5, correction = "cailliez", time_tree = NULL, estimate.all.nodes = FALSE, estimate.tip.values = FALSE, inapplicables.as.missing = FALSE, ancestral.polymorphism.behaviour = "equalp", ancestral.uncertainty.behaviour = "equalp", threshold = 0.01, allow.all.missing = FALSE) {
  
  # Add some top level conditionsl here to check input is valid.
  # Allow other ordination types such as NMDS
  
  # If no tree is supplied:
  if (is.null(time_tree)) {
    
    # Get morphological distances from the cladistic matrix:
    morph_distances <- calculate_morphological_distances(cladistic.matrix, distance.metric = distance.metric, distance.transformation = distance.transformation, polymorphism.behaviour = distance.polymorphism.behaviour, uncertainty.behaviour = distance.uncertainty.behaviour, inapplicable.behaviour = distance.inapplicable.behaviour, character.dependencies = character.dependencies, alpha = alpha)
    
    # Get trimmed distances:
    trimmed.distances <- trim_matrix(morph_distances$DistanceMatrix)

    # If trimming of matrix lead to taxa being removed warn user:
    if (!is.null(trimmed.distances$removed.taxa)) message(paste("The following taxa had to be removed to produce a complete distance matrix:", paste(trimmed.distances$removed.taxa, collapse = ", ")))
    
    # Perform Principal Coordinates Analysis on the data:
    pcoa.results <- ape::pcoa(trimmed.distances$DistMatrix, correction = correction, rn = rownames(trimmed.distances$DistMatrix))
    
  # Case if a tree is included (and a phylomorphospace is requested):
  } else {
      
    # Get ancestral character states:
    ancestral_values <- estimate_ancestral_states(cladistic.matrix = cladistic.matrix, time_tree = time_tree, estimate.all.nodes = estimate.all.nodes, estimate.tip.values = estimate.tip.values, inapplicables.as.missing = inapplicables.as.missing, polymorphism.behaviour = ancestral.polymorphism.behaviour, uncertainty.behaviour = ancestral.uncertainty.behaviour, threshold = threshold, allow.all.missing = allow.all.missing)

    # Get morphological distances from the cladistic matrix:
    morph_distances <- calculate_morphological_distances(ancestral_values, distance.metric = distance.metric, ged.type = ged.type, distance.transformation = distance.transformation, polymorphism.behaviour = distance.polymorphism.behaviour, uncertainty.behaviour = distance.uncertainty.behaviour, inapplicable.behaviour = distance.inapplicable.behaviour, character.dependencies = character.dependencies, alpha = alpha)
    
    # Get trimmed distances:
    trimmed.distances <- trim_matrix(morph_distances$DistanceMatrix, Tree = time_tree)

    # If trimming of matrix lead to taxa or nodes being removed warn user:
    if (!is.null(trimmed.distances$removed.taxa)) message(paste("The following taxa or nodes had to be removed to produce a complete distance matrix:", paste(trimmed.distances$removed.taxa, collapse = ", ")))

    # Store (possibly trimmed) tree ready to be output:
    time_tree <- trimmed.distances$Tree

    # Perform Principal Coordinates Analysis on the data:
    pcoa.results <- ape::pcoa(trimmed.distances$DistMatrix, correction = correction, rn = rownames(trimmed.distances$DistMatrix))

  }
  
  # If corrected vectors exist:
  if (!is.null(pcoa.results$vectors.cor)) {
    
    # Overwrite uncorrected with corrected:
    pcoa.results$vectors <- pcoa.results$vectors.cor
    
    # Remove corrected from output:
    pcoa.results$vectors.cor <- NULL
    
  }
  
  # If corrected trace exists:
  if (!is.null(pcoa.results$trace.cor)) {
    
    # Overwite uncorrected with corrected:
    pcoa.results$trace <- pcoa.results$trace.cor
    
    # Remove corrected from output:
    pcoa.results$trace.cor <- NULL
    
  }
  
  # Compile output:
  output <- c(time_tree = list(time_tree), DistMatrix = list(trimmed.distances$DistMatrix), RemovedTaxa = list(trimmed.distances$RemovedTaxa), pcoa.results)

  # Return output invisibly:
  invisible(output)
  
}
