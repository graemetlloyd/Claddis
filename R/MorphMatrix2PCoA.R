#' Principal Coordinates on a Cladistic Matrix
#'
#' @description
#'
#' Performs Principal Coordinates Analysis (PCoA) on a cladistic matrix.
#'
#' @param CladisticMatrix A vector of mode character representing the tip names for which an ancestor is sought.
#' @param Distance The distance method to use (one of "RED", "GED", "GC", or "MORD" - the default). See \link{MorphDistMatrix} for more details.
#' @param GEDType The type of GED use. Must be one of \code{"Legacy"}, \code{"Hybrid"}, or \code{"Wills"} (the default). See details for an explanation.
#' @param TransformDistances The transformation to apply to distances. See \link{MorphDistMatrix} for details.
#' @param DistPolymorphismBehaviour The distance behaviour for dealing with polymorphisms. Must be one of \code{"mean.difference"}, \code{"min.difference"} (the default), or \code{"random"}. See \link{MorphDistMatrix} for details.
#' @param DistUncertaintyBehaviour The distance behaviour for dealing with uncertainties. Must be one of \code{"mean.difference"}, \code{"min.difference"} (the default), or \code{"random"}. See \link{MorphDistMatrix} for details.
#' @param DistInapplicableBehaviour The behaviour for dealing with inapplicables. Must be one of \code{"missing"} (default), or \code{"HSJ"}. See \link{MorphDistMatrix} for details.
#' @param CharacterDependencies Only relevant if using \code{InapplicableBehaviour = "HSJ"}. Must be a two-column matrix with colnames "DependentCharacter" and "IndependentCharacter" that specifies character hierarchies. See \link{MorphDistMatrix} for details.
#' @param Alpha The alpha value (sensu Hopkins and St John 2018). Only relevant if using \code{InapplicableBehaviour = "HSJ"}. See \link{MorphDistMatrix} for details.
#' @param correction The negative eigenvalue correction to use (one of "lingoes", "none", or "cailliez" - the default). See \link{pcoa} for more details.
#' @param Tree If a phylmorphospace is desired then a tree with root age and branch-lengths must be included.
#' @param EstimateAllNodes If including a tree whether you want to estinate ancestral states for all characters (default is FALSE). See \link{AncStateEstMatrix} for more details.
#' @param EstimateTipValues If including a tree whether you want to estinate missing or polymorphic tip states (default is FALSE). See \link{AncStateEstMatrix} for more details.
#' @param InapplicablesAsMissing See \link{AncStateEstMatrix}.
#' @param AncestralPolymorphismBehaviour Behaviour for dealing with polymorphisms when producing ancestral state estimates - see \link{AncStateEstMatrix}.
#' @param AncestralUncertaintyBehaviour Behaviour for dealing with uncertainties when producing ancestral state estimates - see \link{AncStateEstMatrix}.
#' @param Threshold Threshold for ancestral state estimation of discrete characters - see \link{AncStateEstMatrix} for details.
#'
#' @details
#'
#' Takes a cladistic matrix in the format imported by \link{ReadMorphNexus} and performs Principal Coordinates (Gower 1966) analysis on it.
#'
#' This function is effectively a wrapper for \link{pcoa} from the \link{ape} package and the user is referred there for some of the options (e.g., using the Caillez 1983 approach to avoiding negative eigenvalues).
#'
#' If providing a tree and inferring ancestral states then options to also infer missing or uncertain tips and whether to infer values for all characters at all internal nodes are provided (via \link{AncStateEstMatrix}).
#'
#' Other options within the function concern the distance metric to use and the transformation to be used if selecting a propotional distance (see \link{MorphDistMatrix}).
#'
#' IMPORTANT: The function can remove taxa (or if including a tree, nodes as well) if they lead to an incomplete distance matrix (see \link{TrimMorphDistMatrix}).
#'
#' @return \item{Tree}{The tree (if supplied). Note this may be pruned from the input tree by \link{TrimMorphDistMatrix}.}
#' @return \item{DistMatrix}{The distance matrix. Note this may be pruned by \link{TrimMorphDistMatrix} and thus not include all taxa.}
#' @return \item{RemovedTaxa}{A vector of taxa (or nodes) removed by \link{TrimMorphDistMatrix}. Returns NULL if none are removed.}
#' @return \item{note}{See \link{pcoa}.}
#' @return \item{values}{See \link{pcoa}.}
#' @return \item{vectors}{See \link{pcoa}.}
#' @return \item{trace}{See \link{pcoa}.}
#' @return \item{vectors.cor}{See \link{pcoa}.}
#' @return \item{trace.cor}{See \link{pcoa}.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Cailliez, F., 1983. The analytical solution of the additive constant problem. Psychometrika, 48, 305-308.
#'
#' Gower, J. C., 1966. Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika, 53, 325-338.
#'
#' @examples
#'
#' # Run on Michaux (189) data set with defaults:
#' x <- MorphMatrix2PCoA(Michaux1989)
#'
#' # Show output:
#' x
#'
#' # Generate a (made up) tree:
#' Tree <- rtree(length(rownames(Michaux1989$Matrix_1$Matrix)))
#'
#' # Add taxon names to it:
#' Tree$tip.label <- rownames(Michaux1989$Matrix_1$Matrix)
#'
#' # Set root time by making youngest taxon extant:
#' Tree$root.time <- max(diag(vcv(Tree)))
#'
#' # Run with tree:
#' y <- MorphMatrix2PCoA(Michaux1989, Tree = Tree)
#'
#' # Show new output:
#' y
#'
#' @export MorphMatrix2PCoA
MorphMatrix2PCoA <- function(CladisticMatrix, Distance = "MORD", GEDType = "Wills", TransformDistances = "arcsine_sqrt", DistPolymorphismBehaviour = "min.difference", DistUncertaintyBehaviour = "min.difference", DistInapplicableBehaviour = "missing", CharacterDependencies = NULL, Alpha = 0.5, correction = "cailliez", Tree = NULL, EstimateAllNodes = FALSE, EstimateTipValues = FALSE, InapplicablesAsMissing = FALSE, AncestralPolymorphismBehaviour = "equalp", AncestralUncertaintyBehaviour = "equalp", Threshold = 0.01) {
    
# Add some top level conditionsl here to check input is valid.
  
  # If no tree is supplied:
  if(is.null(Tree)) {
    
    # Get morphological distances from the cladistic matrix:
    morph_distances <- MorphDistMatrix(CladisticMatrix, Distance = Distance, TransformDistances = TransformDistances, PolymorphismBehaviour = DistPolymorphismBehaviour, UncertaintyBehaviour = DistUncertaintyBehaviour, InapplicableBehaviour = DistInapplicableBehaviour, CharacterDependencies = CharacterDependencies, Alpha = Alpha)
    
    # Get trimmed distances:
    trimmed_distances <- TrimMorphDistMatrix(morph_distances$DistanceMatrix)

    # If trimming of matrix lead to taxa being removed warn user:
    if(!is.null(trimmed_distances$removed.taxa)) message(paste("The following taxa had to be removed to produce a complete distance matrix:", paste(trimmed_distances$removed.taxa, collapse = ", ")))
    
    # Perform Principal Coordinates Analysis on the data:
    pcoa_results <- pcoa(trimmed_distances$DistMatrix, correction = correction, rn = rownames(trimmed_distances$DistMatrix))
    
  # Case if a tree is included (and a phylomorphospace is requested):
  } else {
      
    # Get ancestral character states:
    ancestral_values <- AncStateEstMatrix(CladisticMatrix = CladisticMatrix, Tree = Tree, EstimateAllNodes = EstimateAllNodes, EstimateTipValues = EstimateTipValues, InapplicablesAsMissing = InapplicablesAsMissing, PolymorphismBehaviour = AncestralPolymorphismBehaviour, UncertaintyBehaviour = AncestralUncertaintyBehaviour, Threshold = Threshold)

    # Get morphological distances from the cladistic matrix:
    morph_distances <- MorphDistMatrix(ancestral_values, Distance = Distance, GEDType = GEDType, TransformDistances = TransformDistances, PolymorphismBehaviour = DistPolymorphismBehaviour, UncertaintyBehaviour = DistUncertaintyBehaviour, InapplicableBehaviour = DistInapplicableBehaviour, CharacterDependencies = CharacterDependencies, Alpha = Alpha)
    
    # Get trimmed distances:
    trimmed_distances <- TrimMorphDistMatrix(morph_distances$DistanceMatrix, Tree = Tree)

    # If trimming of matrix lead to taxa or nodes being removed warn user:
    if(!is.null(trimmed_distances$removed.taxa)) message(paste("The following taxa or nodes had to be removed to produce a complete distance matrix:", paste(trimmed_distances$removed.taxa, collapse = ", ")))

    # Store (possibly trimmed) tree ready to be output:
    Tree <- trimmed_distances$Tree

    # Perform Principal Coordinates Analysis on the data:
    pcoa_results <- pcoa(trimmed_distances$DistMatrix, correction = correction, rn = rownames(trimmed_distances$DistMatrix))

  }
  
  # Compile output:
  output <- c(list(Tree), list(trimmed_distances$DistMatrix), list(trimmed_distances$RemovedTaxa), pcoa_results)

  # Add variable name for tree:
  names(output)[[1]] <- "Tree"

  # Add variable name for distaance matrix:
  names(output)[[2]] <- "DistMatrix"
  
  # Add variable name for removed taxa and nodes:
  names(output)[[3]] <- "RemovedTaxa"

  # Return output invisibly:
  invisible(output)
  
}
