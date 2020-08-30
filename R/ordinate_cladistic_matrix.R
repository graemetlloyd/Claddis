#' Principal Coordinates on a Cladistic Matrix
#'
#' @description
#'
#' Performs Principal Coordinates Analysis (PCoA) on a cladistic matrix.
#'
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}.
#' @param distance_metric See \link{calculate_morphological_distances}.
#' @param ged_type See \link{calculate_morphological_distances}.
#' @param distance_transformation See \link{calculate_morphological_distances}.
#' @param distance_polymorphism_behaviour See \link{calculate_morphological_distances}.
#' @param distance_uncertainty_behaviour See \link{calculate_morphological_distances}.
#' @param distance_inapplicable_behaviour See \link{calculate_morphological_distances}.
#' @param character_dependencies See \link{calculate_morphological_distances}.
#' @param alpha See \link{calculate_morphological_distances}.
#' @param correction The negative eigenvalue correction to use (one of "lingoes", "none", or "cailliez" - the default). See \link{pcoa} for more details.
#' @param time_tree If a phylmorphospace is desired then a tree with root age and branch-lengths must be included.
#' @param estimate_all_nodes See \link{estimate_ancestral_states}.
#' @param estimate_tip_values See \link{estimate_ancestral_states}.
#' @param inapplicables_as_missing See \link{estimate_ancestral_states}.
#' @param ancestral_polymorphism_behaviour See \link{estimate_ancestral_states}.
#' @param ancestral_uncertainty_behaviour See \link{estimate_ancestral_states}.
#' @param threshold See \link{estimate_ancestral_states}.
#' @param all_missing_allowed See \link{estimate_ancestral_states}.
#'
#' @details
#'
#' Takes a cladistic matrix in the format imported by \link{read_nexus_matrix} and performs Principal Coordinates (Gower 1966) analysis on it.
#'
#' This function is effectively a wrapper for the pipeline:
#'
#' \link{estimate_ancestral_states} -> \link{calculate_morphological_distances} -> \link{pcoa}
#'
#' With the first part being optional (if wanting a phylomorphospace) and the latter coming from the \link{ape} package (the user is referred there for some of the options, e.g., using the Caillez 1983 approach to avoiding negative eigenvalues). (See Lloyd 2016 for more on disparity pipelines.)
#'
#' If providing a tree and inferring ancestral states then options to also infer missing or uncertain tips and whether to infer values for all characters at all internal nodes are provided by the \link{estimate_ancestral_states} part.
#'
#' Other options within the function concern the distance metric to use and the transformation to be used if selecting a propotional distance (see \link{calculate_morphological_distances}).
#'
#' IMPORTANT: The function can remove taxa (or if including a tree, nodes as well) if they lead to an incomplete distance matrix (see \link{trim_matrix} for more details).
#'
#' @return
#'
#' \item{time_tree}{The tree (if supplied). Note this may be pruned from the input tree by \link{trim_matrix}.}
#' \item{distance_matrix}{The distance matrix. Note this may be pruned by \link{trim_matrix} and thus not include all taxa.}
#' \item{removed_taxa}{A vector of taxa and/or nodes removed by \link{trim_matrix}. Returns NULL if none were removed.}
#' \item{note}{See \link{pcoa}.}
#' \item{values}{See \link{pcoa}.}
#' \item{vectors}{See \link{pcoa}. Note: this will be the same as \code{vectors.cor} from the \link{pcoa} output if a correction was applied.}
#' \item{trace}{See \link{pcoa}. Note: this will be the same as \code{trace.cor} from the \link{pcoa} output if a correction was applied.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{assign_taxa_to_bins}, \link{plot_chronophylomorphospace}, \link{plot_morphospace_stack}, \link{plot_morphospace}, \link{plot_multi_morphospace}
#'
#' @references
#'
#' Cailliez, F., 1983. The analytical solution of the additive constant problem. \emph{Psychometrika}, \bold{48}, 305-308.
#'
#' Gower, J. C., 1966. Some distance properties of latent root and vector methods used in multivariate analysis. \emph{Biometrika}, \bold{53}, 325-338.
#'
#' @examples
#'
#' # Run on Michaux (1989) data set with default settings:
#' x <- ordinate_cladistic_matrix(cladistic_matrix = michaux_1989)
#'
#' # Show entire output:
#' x
#'
#' # Generate a (made up) tree:
#' time_tree <- ape::rtree(n = length(x = rownames(x = michaux_1989$matrix_1$matrix)))
#'
#' # Add taxon names to it:
#' time_tree$tip.label <- rownames(x = michaux_1989$matrix_1$matrix)
#'
#' # Set root time by making youngest taxon extant:
#' time_tree$root.time <- max(diag(x = ape::vcv(phy = time_tree)))
#'
#' # Run with tree:
#' y <- ordinate_cladistic_matrix(cladistic_matrix = michaux_1989, time_tree = time_tree)
#'
#' # Show new output:
#' y
#' @export ordinate_cladistic_matrix
ordinate_cladistic_matrix <- function(cladistic_matrix, distance_metric = "mord", ged_type = "wills", distance_transformation = "arcsine_sqrt", distance_polymorphism_behaviour = "min_difference", distance_uncertainty_behaviour = "min_difference", distance_inapplicable_behaviour = "missing", character_dependencies = NULL, alpha = 0.5, correction = "cailliez", time_tree = NULL, estimate_all_nodes = FALSE, estimate_tip_values = FALSE, inapplicables_as_missing = FALSE, ancestral_polymorphism_behaviour = "equalp", ancestral_uncertainty_behaviour = "equalp", threshold = 0.01, all_missing_allowed = FALSE) {

  # Add some top level conditionsl here to check input is valid.
  # Allow other ordination types such as NMDS.
  # Allow the "bad" phylomorphospace type (ordination axis continuous ancestors).
  # Add scree values to output and use in plotting functions.

  # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
  if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")

  # If no tree is supplied:
  if (is.null(time_tree)) {

    # Get morphological distances from the cladistic matrix:
    morphological_distances <- calculate_morphological_distances(cladistic_matrix = cladistic_matrix, distance_metric = distance_metric, distance_transformation = distance_transformation, polymorphism_behaviour = distance_polymorphism_behaviour, uncertainty_behaviour = distance_uncertainty_behaviour, inapplicable_behaviour = distance_inapplicable_behaviour, character_dependencies = character_dependencies, alpha = alpha)

    # Get trimmed distances:
    trimmed_distances <- trim_matrix(morphological_distances$distance_matrix)

    # If trimming of matrix lead to taxa being removed warn user:
    if (!is.null(trimmed_distances$removed_taxa)) message(paste("The following taxa had to be removed to produce a complete distance matrix:", paste(trimmed_distances$removed_taxa, collapse = ", ")))

    # Perform Principal Coordinates Analysis on the data:
    pcoa_results <- ape::pcoa(trimmed_distances$distance_matrix, correction = correction, rn = rownames(x = trimmed_distances$distance_matrix))

    # Case if a tree is included (and a phylomorphospace is requested):
  } else {

    # Get ancestral character states:
    ancestral_values <- estimate_ancestral_states(cladistic_matrix = cladistic_matrix, time_tree = time_tree, estimate_all_nodes = estimate_all_nodes, estimate_tip_values = estimate_tip_values, inapplicables_as_missing = inapplicables_as_missing, polymorphism_behaviour = ancestral_polymorphism_behaviour, uncertainty_behaviour = ancestral_uncertainty_behaviour, threshold = threshold, all_missing_allowed = all_missing_allowed)

    # Get morphological distances from the cladistic matrix:
    morphological_distances <- calculate_morphological_distances(cladistic_matrix = ancestral_values, distance_metric = distance_metric, ged_type = ged_type, distance_transformation = distance_transformation, polymorphism_behaviour = distance_polymorphism_behaviour, uncertainty_behaviour = distance_uncertainty_behaviour, inapplicable_behaviour = distance_inapplicable_behaviour, character_dependencies = character_dependencies, alpha = alpha)

    # Get trimmed distances:
    trimmed_distances <- trim_matrix(morphological_distances$distance_matrix, tree = time_tree)

    # If trimming of matrix lead to taxa or nodes being removed warn user:
    if (!is.null(trimmed_distances$removed_taxa)) message(paste("The following taxa or nodes had to be removed to produce a complete distance matrix:", paste(trimmed_distances$removed_taxa, collapse = ", ")))

    # Store (possibly trimmed) tree ready to be output:
    time_tree <- trimmed_distances$tree

    # Perform Principal Coordinates Analysis on the data:
    pcoa_results <- ape::pcoa(trimmed_distances$distance_matrix, correction = correction, rn = rownames(x = trimmed_distances$distance_matrix))
  }

  # If corrected vectors exist:
  if (!is.null(pcoa_results$vectors.cor)) {

    # Overwrite uncorrected with corrected:
    pcoa_results$vectors <- pcoa_results$vectors.cor

    # Remove corrected from output:
    pcoa_results$vectors.cor <- NULL
  }

  # If corrected trace exists:
  if (!is.null(pcoa_results$trace.cor)) {

    # Overwite uncorrected with corrected:
    pcoa_results$trace <- pcoa_results$trace.cor

    # Remove corrected from output:
    pcoa_results$trace.cor <- NULL
  }

  # Return compiled output:
  invisible(c(time_tree = list(time_tree), distance_matrix = list(trimmed_distances$distance_matrix), removed_taxa = list(trimmed_distances$removed_taxa), pcoa_results))
}
